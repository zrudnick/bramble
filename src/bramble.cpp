// Bramble v1.0.0
// Main file

#include "bramble.h"
#include "bam.h"
#include "reads.h"
#include "tree.h"
#include "proc_mem.h" 			// for GMEMTRACE
#include "GThreads.h" 			// for THREADING_ENABLED
#include "htslib/sam.h"
#include "CLI/CLI11.hpp"
#include <numeric>

#define VERSION "1.0.0"

#define USAGE "Bramble v" VERSION " usage:\n\n\
bramble <in.bam ..> [-G <guide_gff>] [-o <out.gtf>] [-p <cpus>]\n\
 [-v] [-g <bdist>] [-u] [-L] [-h] [--rc] [--fw]\n\
 \n\
Assemble RNA-Seq alignments into potential transcripts.\n\
Options:\n\
 --version : print just the version at stdout and exit\n\
 --fr : assume stranded library fw-firststrand\n\
 --rf : assume stranded library fw-secondstrand\n\
 -G reference annotation to use for guiding the assembly process (GTF/GFF)\n\
 -o output path/file name for the assembled transcripts GTF (default: stdout)\n\
 -L BAM file contains long reads \n\
 -s minimum reads per bp coverage to consider for single-exon transcript\n\
    (default: 4.75)\n\
 -v verbose (log bundle processing details)\n\
 -g maximum gap allowed between read mappings (default: 50)\n\
 -p number of threads (CPUs) to use (default: 1)\n\
 -u no multi-mapping correction (default: correction enabled)\n\
 -h print this usage message and exit\n\
"

FILE* f_out=NULL;
FILE* c_out=NULL;
FILE* dbg_out=NULL;

char* bam_file_in;
const char* bam_file_out;

GStr outfname;
GStr out_dir;
GStr tmp_fname;						// Output file name (from command line, -o)
GStr guide_gff; 					// Reference annotation, -G
GFastaDb* gfasta = NULL;

bool debugMode = false;
bool verbose = false;
bool longreads = false;				// Long reads were used in mapping
bool use_fasta = false;				// Use FASTA for more precise transcriptome mapping
bool use_sc = false;				// Single-cell data was used in mapping
bool fr_strand = false;				// Read 1 is on forward strand, Read 2 is on reverse strand
bool rf_strand = false;				// Read 1 is on reverse strand, Read 2 is on forward strand
int num_cpus = 1;					// Threads, -p

uint bundledist = 50;  				// Reads at what distance should be considered part of separate bundles
uint runoffdist = 200;
uint junctionsupport = 100; 	// anchor length for junction to be considered well supported <- consider shorter??

GffNames* guide_seq_names = NULL; 	// Used as a dictionary for reference sequence names and ids
int num_ref_seqs = 0; 				// Number of reference sequences found in the guides file

// For GMEMTRACE
double maxMemRS=0;
double maxMemVM=0;
GStr maxMemBundle;

#ifndef NOTHREADS
// Threading: single producer, multiple consumers
// Main thread/program is always loading the producer

GMutex data_mutex; 					// Manage availability of data records ready to be loaded by main thread
GVec<int> clear_data_pool; 			// Indexes of data bundles cleared for loading by main thread (clear data pool)
GConditionVar have_bundles; 		// Will notify a thread that a bundle was loaded in the ready queue
                           			// (or that no more bundles are coming)
int bundle_work = 1; 				// Bit 0 set if bundles are still being prepared (BAM file not exhausted yet)
                  					// Bit 1 set if there are Bundles ready in the queue

GMutex wait_mutex; 					// Controls threads_waiting (idle threads counter)
int threads_waiting; 				// Idle worker threads
GConditionVar haveThreads; 			// Will notify the bundle loader when a thread
                          			// Is available to process the currently loaded bundle

GConditionVar have_clear; 			// Will notify when bundle buf space available
GMutex queue_mutex; 				// Controls bundle_queue and bundles access
GFastMutex printMutex; 				// For writing the output to file
GFastMutex logMutex; 				// Only when verbose - to avoid mangling the log output
GFastMutex bam_reading_mutex;
GFastMutex printCovMutex;
#endif

bool no_more_bundles = false;
bool has_more_bundles(); 				// Thread-safe retrieves no_more_bundles
void wait_for_bundles(); 				// Sets no_more_bundles to true

void process_options(GArgs& args);
char* sprint_time();
void process_bundle(BundleData* bundle, BamIO* io);
bool no_threads_waiting();
void worker_thread(GThreadData& td); 	    // Thread function
int wait_for_data(BundleData* bundles);		//prepare the next free bundle for loading

void printTime(FILE* f) {
	time_t ltime; /* calendar time */
	ltime=time(NULL);
	struct tm *t=localtime(&ltime);
	fprintf(f, "[%02d/%02d %02d:%02d:%02d]",t->tm_mon+1, t->tm_mday,
			t->tm_hour, t->tm_min, t->tm_sec);
}

// Print the time to a string buffer
char* sprint_time() {
	static char sbuf[32];
	time_t ltime; /* calendar time */
	ltime=time(NULL);
	struct tm *t=localtime(&ltime);
	sprintf(sbuf, "%02d_%02d_%02d:%02d:%02d",t->tm_mon+1, t->tm_mday,
			t->tm_hour, t->tm_min, t->tm_sec);
	return(sbuf);
}

// Process input parameters
void process_options(GArgs& args) {
	// Help
	if (args.getOpt('h') || args.getOpt("help")) {
		fprintf(stdout,"%s",USAGE);
		exit(0);
	}

	// Version
	if (args.getOpt("version")) {
	fprintf(stdout,"%s\n",VERSION);
	exit(0);
	}

	// FR (Forward-Reverse)
	// First read is on forward strand, second read is on reverse strand
	if (args.getOpt("fr")) fr_strand=true;

	// RF (Reverse-Forward)
	// First read is on reverse strand, second read is on forward strand
	if (args.getOpt("rf")) {
		rf_strand=true;
		if(fr_strand) GError("Error: --fr and --rf options are incompatible.\n");
	}

	// Debug Mode
	debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);

	// Verbose
	verbose=(args.getOpt('v')!=NULL);
	if (verbose) {
		fprintf(stderr, "Running StringTie " VERSION ". Command line:\n");
		args.printCmdLine(stderr);
	}

	GStr s;

	// Number of Threads, -p
	s=args.getOpt('p');
	if (!s.is_empty()) {
		num_cpus=s.asInt();
		if (num_cpus<=0) num_cpus=1;
	}

	// Bundle Distance, -g
	s=args.getOpt('g');
	if (!s.is_empty()) {
		bundledist=s.asInt();
		if (bundledist>runoffdist) runoffdist = bundledist;
	}

	// Reference Annotation
	if (args.getOpt('G')) {
	guide_gff=args.getOpt('G');
	if (!fileExists(guide_gff.chars())>1) 
	    GError("Error: reference annotation file (%s) not found.\n",
				guide_gff.chars());
	}

	s = args.getOpt('S');
	if (!s.is_empty()) {
		gfasta = new GFastaDb(s.chars());
		use_fasta = true;
	}

    s=args.getOpt('C');
    if (!s.is_empty()) {
        c_out=fopen(s.chars(), "w");
        if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
    }

	// Long Reads
	longreads=(args.getOpt('L')!=NULL);
	if (longreads) bundledist = 0;

	int numbam=args.startNonOpt();

	if (numbam < 1) {
		GMessage("%s\nError: no input file provided!\n",USAGE);
		exit(1);
	}

	if (guide_gff == NULL) {
		GMessage("%s\nError: no input file provided!\n",USAGE);
		exit(1);
	}

	char* ifn=NULL;
	int i = 0;
	while ( (ifn=args.nextNonOpt())!=NULL) {
		//input alignment files
		if (i == 0) bam_file_in = ifn;
	}

	// Create output path

	// Output path, -o
	tmp_fname = args.getOpt('o');
	bam_file_out = tmp_fname;
	GMessage("Output file: %s\n", bam_file_out);

	if (bam_file_out == NULL) {
		GMessage("%s\nError: no output file name provided!\n",USAGE);
		exit(1);
	}
	
	// outfname="stdout";
	// out_dir="./";
	// if (!tmp_fname.is_empty() && tmp_fname!="-") {
	// 	if (tmp_fname[0]=='.' && tmp_fname[1]=='/')
	// 		tmp_fname.cut(0,2);
	// 	outfname=tmp_fname;
	// 	int pidx=outfname.rindex('/');
	// 	if (pidx>=0) {//path given
	// 		out_dir=outfname.substr(0,pidx+1);
	// 		tmp_fname=outfname.substr(pidx+1);
	// 	}
	// }
	// else { // stdout
	// 	tmp_fname=outfname;
	// 	char *stime=sprint_time();
	// 	tmp_fname.tr(":","-");
	// 	tmp_fname+='.';
	// 	tmp_fname+=stime;
	// }
	// if (out_dir!="./") {
	// 	if (fileExists(out_dir.chars())==0) {
	// 		//directory does not exist, create it
	// 		if (Gmkdir(out_dir.chars()) && !fileExists(out_dir.chars())) {
	// 			GError("Error: cannot create directory %s!\n", out_dir.chars());
	// 		}
	// 	}
	// }
 
	// if(outfname != "stdout") {
	// 	GStr bam_file_out = outfname;
	// 	if (bam_file_out==NULL) GError("Error creating output file %s\n", bam_file_out.chars());
	// }
}

// Check if there are more bundles
bool has_more_bundles() {
	if (THREADING_ENABLED) GLockGuard<GFastMutex> lock(bam_reading_mutex);
	return !no_more_bundles;
}

// Wait for threads and set bundle status
void wait_for_bundles() {

	if (THREADING_ENABLED) {

		// BAM Reading Mutex
		bam_reading_mutex.lock();
		no_more_bundles = true;
		bam_reading_mutex.unlock();

		// Queue Mutex
		queue_mutex.lock();
		bundle_work &= ~(int)0x01; //clear bit 0;
		queue_mutex.unlock();

		bool are_threads_waiting=true;
		while (are_threads_waiting) {

			// Wait Mutex
			wait_mutex.lock();
			are_threads_waiting = (threads_waiting > 0);
			wait_mutex.unlock();

			if (are_threads_waiting) {

				have_bundles.notify_all();
				current_thread::sleep_for(1);

				// Wait Mutex
				wait_mutex.lock();
				are_threads_waiting = (threads_waiting > 0);
				wait_mutex.unlock();

				current_thread::sleep_for(1);
			}
		}
	} 
	
	else no_more_bundles=true;
}

// Process current bundle
void process_bundle(BundleData* bundle, BamIO* io) {

	if (verbose) {
		if (THREADING_ENABLED) GLockGuard<GFastMutex> lock(logMutex);

		printTime(stderr);
		GMessage(">bundle %s:%d-%d [%lu alignments (%d distinct), %d junctions, %d guides] begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->reads.Count(), bundle->junction.Count(),
				bundle->guides.Count());
		if (GMEMTRACE) {
			double vm;
			double rsm;
			get_mem_usage(vm, rsm);
			GMessage("\t\tstart memory usage: %6.1fMB\n",rsm/1024);
			if (rsm>maxMemRS) {
				maxMemRS=rsm;
				maxMemVM=vm;
				maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->reads.Count());
			}
		}
	}

	// This is where reads are added to new BAM file
	convert_reads(bundle, io);

	if (verbose) {
		if (THREADING_ENABLED) GLockGuard<GFastMutex> lock(logMutex);
		
		printTime(stderr);
		GMessage("^bundle %s:%d-%d done (%d processed potential transcripts).\n",bundle->refseq.chars(), bundle->start, bundle->end);
		if (GMEMTRACE) {
			double vm;
			double rsm;
			get_mem_usage(vm, rsm);
			GMessage("\t\tfinal memory usage: %6.1fMB\n",rsm/1024);
			if (rsm>maxMemRS) {
				maxMemRS=rsm;
				maxMemVM=vm;
				maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->reads.Count());
			}
		}
	}
	bundle->Clear();
}

// Check that there aren't any threads waiting
bool no_threads_waiting() {
	wait_mutex.lock();
	int threads = threads_waiting;
	wait_mutex.unlock();
	return (threads < 1);
}

// Worker thread waits for incoming bundle, calls process_bundle
void worker_thread(GThreadData& td) {
	WorkerArgs* args = (WorkerArgs*)(td.udata);
	GPVec<BundleData>* bundle_queue = args->bundle_queue;
	BamIO* io = args->io;

	// Wait for a ready bundle in the queue, until there is no hope for incoming bundles
	queue_mutex.lock(); //enter wait-for-notification loop

	while (bundle_work) {
		wait_mutex.lock();
		threads_waiting++;
		queue_mutex.unlock();
		wait_mutex.unlock();
		haveThreads.notify_one(); //in case main thread is waiting
		current_thread::yield();
		queue_mutex.lock();
		while (bundle_work && bundle_queue->Count()==0) {
			have_bundles.wait(queue_mutex); //unlocks queue_mutex and wait until notified
					//when notified, locks queue_mutex and resume
		}

		wait_mutex.lock();
		if (threads_waiting > 0) threads_waiting--;
		wait_mutex.unlock();

		BundleData* readyBundle = NULL;
		if ((bundle_work & 0x02)!=0) { //is bit 1 set?
			readyBundle = bundle_queue->Pop();

			// Make sure bundle is not NULL
			if (readyBundle != NULL) {
				if (bundle_queue->Count()==0)
				bundle_work &= ~(int)0x02; //clear bit 1 (queue is empty)
			
				// Queue Mutex
				queue_mutex.unlock();
				process_bundle(readyBundle, io);

				// Data Mutex
				data_mutex.lock();
				clear_data_pool.Push(readyBundle->idx);
				data_mutex.unlock();
				
				have_clear.notify_one(); //inform main thread
				current_thread::yield();
				
				queue_mutex.lock();
			}
		}
	}
	queue_mutex.unlock();
}

// Prepare the next available bundle slot for loading
int wait_for_data(BundleData* bundles) {
	int bidx=-1;

	data_mutex.lock();
	while (clear_data_pool.Count()==0) {
		have_clear.wait(data_mutex);
	}
	bidx=clear_data_pool.Pop();
	if (bidx>=0) bundles[bidx].status = BundleStatus::BUNDLE_STATUS_LOADING;
	data_mutex.unlock();

	return bidx;
}

void rc_updateExonCounts(const RC_ExonOvl& exonovl, int nh) {
  	//this only gets read overlaps > 5bp and otherwise filtered in evalReadAln()
  	exonovl.feature->rcount++;
  	if (nh>1) {
		exonovl.feature->mrcount += (1.0/nh);
	  	exonovl.feature->movlcount += ((double)exonovl.ovlen/nh);
  	}
  	else { // nh<=1
	  	exonovl.feature->mrcount++;
	  	exonovl.feature->movlcount += exonovl.ovlen;
	  	exonovl.feature->ucount++;
 	 }
}

bool BundleData::evalReadAln(GReadAlnData& alndata, char& xstrand) {
    if (rc_data == NULL) return false; // no guides available for this read's region

    GSamRecord* brec = alndata.brec;
    int mate_pos = brec->mate_start();
    int nh = alndata.nh;
    if ((int)brec->end < rc_data->lmin || (int)brec->start > rc_data->rmax) return false; //hit outside coverage area
    if (rc_data->g_exons.Count() == 0 || rc_data->g_tdata.Count() == 0) return false; //nothing to do without transcripts

    // Check this read alignment against guide exons and introns
    char strandbits = 0;
    bool overlaps_guide = false;
    bool is_in_guide = true; // exons and junctions are in reference transcripts but they might be in different guides
    
    for (int i = 0; i < brec->exons.Count(); i++) {
        GArray<RC_ExonOvl> exonOverlaps(true, true); //overlaps sorted by decreasing length

        if (rc_data->findOvlExons(exonOverlaps, brec->exons[i].start, brec->exons[i].end, xstrand, mate_pos)) {
            overlaps_guide = true;
            int max_ovl = exonOverlaps[0].ovlen;
             
            if (is_in_guide && (uint)max_ovl < brec->exons[i].len()) is_in_guide = false;

            for (int k = 0; k < exonOverlaps.Count(); ++k) {
                //if (exonOverlaps[k].ovlen < 5) break; //ignore very short overlaps
                if (k && (exonOverlaps[k].mate_ovl < exonOverlaps[0].mate_ovl || exonOverlaps[k].ovlen + 5 < max_ovl))
                    break; //ignore further overlaps after a mate matched or if they are shorter than max_overlap-5

                if (exonOverlaps[k].feature->strand == '+') strandbits |= 0x01;
                else if (exonOverlaps[k].feature->strand == '-') strandbits |= 0x02;

                rc_updateExonCounts(exonOverlaps[k], nh);
            }
        }

        // Intron processing
        if (i > 0) {
            int j_l = brec->exons[i-1].end + 1;
            int j_r = brec->exons[i].start - 1;
            RC_Feature* ri = rc_data->findIntron(j_l, j_r, xstrand);
            alndata.juncs.Add(new CJunction(j_l, j_r)); //don't set strand, etc. for now
            if (ri) { //update guide intron counts
                ri->rcount++;
                ri->mrcount += (nh > 1) ? (1.0/nh) : 1;
                if (nh == 1) ri->ucount++;
                alndata.juncs.Last()->guide_match = 1;
            }
            else is_in_guide = false;
        }
    }

    if (xstrand == '.' && strandbits && strandbits < 3) {
        xstrand = (strandbits == 1) ? '+' : '-';
    }

    return overlaps_guide;
}

int gseqstat_cmpName(const pointer p1, const pointer p2) {
	return strcmp(((GSeqStat*)p1)->gseqname, ((GSeqStat*)p2)->gseqname);
}

// Main function
int main(int argc, char* argv[]) {


 /*
 --fr : assume stranded library fw-firststrand\n\
 --rf : assume stranded library fw-secondstrand\n\
 -G reference annotation to use for guiding the assembly process (GTF/GFF)\n\
 -o output path/file name for the assembled transcripts GTF (default: stdout)\n\
 -L BAM file contains long reads \n\
 -s minimum reads per bp coverage to consider for single-exon transcript\n\
    (default: 4.75)\n\
 -v verbose (log bundle processing details)\n\
 -g maximum gap allowed between read mappings (default: 50)\n\
 -p number of threads (CPUs) to use (default: 1)\n\
 -u no multi-mapping correction (default: correction enabled)\n\
 -h print this usage message and exit\n\
 */


	// Process arguments
	/*
	GArgs args(argc, argv,
   	"debug;help;version;mix;ref=;cram-ref=cds=;keeptmp;rseq=;ptf=;bam;fr;rf;merge;"
   	"exclude=zihvteuLRx:n:j:s:D:G:C:S:l:m:o:a:j:c:f:p:g:P:M:Bb:A:E:F:T:");
 	args.printError(USAGE, true);
 	process_options(args);
 	*/

	
	CLI::App app{"A program to project spliced genomic alignments onto the transcriptome", "Bramble"};
	std::string gff;
	std::string bam_file;
	std::string input_bam;
	app.add_option("in.bam", input_bam, "input bam file")->required();
 	app.add_flag("--fr", fr_strand, "assume stranded library fw-firststrand");
 	app.add_flag("--rf", rf_strand, "assume stranded library fw-secondstrand");
 	app.add_flag("-v", verbose, "verbose (log processing details)");
 	app.add_option("-G", gff, "reference annotation to use for guiding the assembly process (GTF/GFF)")->check(CLI::ExistingFile);
 	app.add_option("-o", bam_file, "output path/file name for the projected alignments")->required();
	app.add_option("-p", num_cpus, "number of threads (CPUs) to use")->default_val(1);

	CLI11_PARSE(app, argc, argv);

	bam_file_in = new char[input_bam.size()+1];
	std::strcpy(bam_file_in, input_bam.c_str());
	guide_gff = gff.c_str();
	bam_file_out = bam_file.c_str();

	GVec<GRefData> refguides;                   		// Plain vector with transcripts for each chromosome

	// Table indexes for raw counts data
	GPVec<RC_TData> guides_RC_tdata(true);     	 		// Raw count data for all guide transcripts
	GPVec<RC_Feature> guides_RC_exons(true);    		// Raw count data for all guide exons
	GPVec<RC_Feature> guides_RC_introns(true);  		// Raw count data for all guide introns

	// Setup debug options
	#ifdef DEBUGPRINT
	verbose = true;
	#endif

	const char* ERR_BAM_SORT = "\nError: the input alignment file is not sorted!\n";

	// Create SAM header file
	const GStr header_path = "tx_header.sam";
    FILE* header_file = fopen(header_path.chars(), "w");
    if (header_file == NULL) GError("Error creating file: %s\n", header_path.chars());
    fprintf(header_file, "@HD\tVN:1.0\tSO:coordinate\n");
	
	if (verbose) {
		printTime(stderr);
		GMessage(" Loading reference annotation (guides)..\n");
	}
	
	// Open GFF/GTF file
	FILE* f = fopen(guide_gff.chars(), "r");
	if (f == NULL) GError("Error: could not open reference annotation file (%s)!\n", guide_gff.chars());

	// GffReader: transcripts only, sort by location
	GffReader gffreader(f, true, true);           // Loading only recognizable transcript features
	gffreader.setRefAlphaSorted();                // Alphabetical sorting of RefSeq IDs
	gffreader.showWarnings(verbose);

	// keep attributes, merge close exons, no exon attributes
	// merge_close_exons must be false for correct transcriptome header construction!
	gffreader.readAll(true, false, false);
	num_ref_seqs = gffreader.gseqtable.Count();
	if (num_ref_seqs == 0 || gffreader.gflst.Count() == 0) {
		GError("Error: could not any valid reference transcripts in %s (invalid GTF/GFF file?)\n", guide_gff.chars());
	}
	//gffreader.gseqStats.Sort(gseqstat_cmpName);
	refguides.setCount(num_ref_seqs); 	//maximum gseq_id

	// Process each transcript for guide loading and header generation
	uint curr_tid = 0;
	int last_refid = -1;

	for (int i = 0; i < gffreader.gflst.Count(); i++) {
		GffObj* guide = gffreader.gflst[i];

		// Chromosome switch
		if (last_refid != guide->gseq_id) {
		   last_refid = guide->gseq_id;
	   }

		// Sanity check: make sure there are no exonless "genes" or other
		if (guide->exons.Count()==0) {
			if (verbose)
				GMessage("Warning: exonless GFF %s feature with ID %s found, added implicit exon %d-%d.\n",
							guide->getFeatureName(), guide->getID(), guide->start, guide->end);
			guide->addExon(guide->start, guide->end); //should never happen!
		}

		// Header generation: calculate spliced transcript length by summing exon lengths
        int transcript_len = 0;
        for (int j = 0; j < guide->exons.Count(); j++) {
            if (guide->exons[j] != NULL) {
                transcript_len += guide->exons[j]->end - guide->exons[j]->start + 1;
            }
        }
        
        // Write to SAM header if we have valid transcript info
        auto id = guide->getID();
        if (id != NULL && transcript_len > 0) {
            fprintf(header_file, "@SQ\tSN:%s\tLN:%d\n", id, transcript_len);
        }

		// Always keep a RC_TData pointer around, with additional info about guides
		RC_TData* tdata = new RC_TData(*guide, ++curr_tid);
		guide->uptr = tdata;
		guides_RC_tdata.Add(tdata);
		GRefData& grefdata = refguides[guide->gseq_id];
	   	grefdata.add(&gffreader, guide); //transcripts already sorted by location

	}

	fprintf(header_file, "@CO\tGenerated from GTF: %s\n", guide_gff.chars());
    fclose(header_file);

	if (verbose) {
		printTime(stderr);
		GMessage("Transcriptome header created\n");
		GMessage(" %d reference transcripts loaded\n", gffreader.gflst.Count());
	}

	// Start BAM IO
	BamIO* io = new BamIO(bam_file_in, bam_file_out, header_path);
	io->start();

	guide_seq_names = GffObj::names; 		// Might have been populated already by gff data
	gffnames_ref(guide_seq_names);  		// Initialize the names collection if not guided

	// Process input BAM file
	// Bundle reads into overlapping groups

	GHash<int> hashread;      			// read_name:pos:hit_index => readlist index		used in mate pairing?
	GList<GffObj>* guides = NULL; 		// List of transcripts on a specific reference
	uint curr_bundle_start = 0;
	uint curr_bundle_end = 0;
	int first_possible_overlap = 0;
	int last_guide_idx = -1;
	int total_guides = 0;

	GStr last_ref_name;
	int last_ref_id = -1; 			//last seen ref_id

#ifndef NOTHREADS
	GMessage("THREADING_ENABLED is on\n");
	
	#define DEF_TSTACK_SIZE 8388608
	size_t def_stack_size = DEF_TSTACK_SIZE;

	#ifdef _GTHREADS_POSIX_
	size_t tstack_size = GThread::defaultStackSize();
	if (tstack_size < DEF_TSTACK_SIZE) def_stack_size = DEF_TSTACK_SIZE;
	if (verbose) {
		if (tstack_size < def_stack_size) GMessage("Default stack size for threads: %d (increased to %d)\n", tstack_size, def_stack_size);
		else GMessage("Default stack size for threads: %d\n", tstack_size);
	}
	#endif

	GThread* threads = new GThread[num_cpus]; 								// Threads for processing bundles
	GPVec<BundleData>* bundle_queue = new GPVec<BundleData>(false); 		// Queue for holding loaded bundles
	BundleData* bundles = new BundleData[num_cpus + 1]; 					// redef with more bundles

	clear_data_pool.setCapacity(num_cpus + 1);

	WorkerArgs* worker_args = new WorkerArgs(bundle_queue, io);

	// Use CPUs to run worker threads
	// Calls worker_thread -> process_bundle -> infer_transcripts
	for (int b = 0; b < num_cpus; b++) {
		threads[b].kickStart(worker_thread, (void*) worker_args, def_stack_size);
		bundles[b+1].idx = b + 1;
		clear_data_pool.Push(b);
	}
	BundleData* bundle = &(bundles[num_cpus]);

	// Otherwise just put everything into the same bundle
#else
	BundleData bundles[1];
	BundleData* bundle = &(bundles[0]);
#endif
	
	GSamRecord* brec = nullptr;
	bool more_alignments = true;
	int prev_pos = 0;

	while (more_alignments) {
		bool new_chr = false;
		int read_start_pos = 0;
		const char* ref_name = NULL;
		char splice_strand = 0;
		int nh = 1;
		int hi = 0;
		int ref_id = last_ref_id;  //current chr id
		bool new_bundle = false;
		
		if ((brec = io->next()) != NULL) {
			
			if (brec->isUnmapped()) continue;
			if (brec->start < 1 || brec->mapped_len < 10) {
				if (verbose) GMessage("Warning: invalid mapping found for read %s (position=%d, mapped length=%d)\n",
						brec->name(), brec->start, brec->mapped_len);
				continue;
			}

			ref_name = brec->refName();

			// Determine splice strand
			splice_strand = brec->spliceStrand(); // tagged strand gets priority
			
			if ((splice_strand == '.') && (fr_strand || rf_strand)) { // set strand if stranded library

				// Read is paired
				if (brec->isPaired()) {
					if (brec->pairOrder() == 1) { // first read in pair
						if ((rf_strand && brec->revStrand()) || (fr_strand && !brec->revStrand())) 
							 splice_strand = '+';
						else splice_strand = '-';
					}
					else {
						if ((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) 
							 splice_strand = '-';
						else splice_strand = '+';
					}

				// Read isn't paired
				} else {
					if ((rf_strand && brec->revStrand()) || (fr_strand && !brec->revStrand())) 
						 splice_strand = '+';
					else splice_strand = '-';
				}

			}

			if (ref_name == NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");
			read_start_pos = brec->start; //BAM is 0 based, but GBamRecord makes it 1-based

			// Check if ref name is new (a different chromosome)
			new_chr = (last_ref_name.is_empty() || last_ref_name != ref_name); 	// chromosome has changed
			if (new_chr) {
				ref_id = guide_seq_names->gseqs.addName(ref_name);
				
				if (ref_id >= num_ref_seqs) {
					if (verbose) {
						GMessage("WARNING: no reference transcripts found for genomic sequence \"%s\"! (mismatched reference names?)\n",
								ref_name); }
				}
				
				prev_pos = 0;
			}
			if (read_start_pos < prev_pos) GError("%s\nread %s (start %d) found at position %d on %s when prev_pos=%d\n",
				ERR_BAM_SORT, brec->name(), brec->start, read_start_pos, ref_name, prev_pos);
			prev_pos = read_start_pos;

			nh = brec->tag_int("NH");		// Number of reported alignments that contain the query in the current record
			if (nh == 0) nh = 1;
			hi = brec->tag_int("HI");		// Query hit index

			if ((!new_chr) && (curr_bundle_end > 0) && (read_start_pos > (curr_bundle_end + (int)runoffdist))) new_bundle = true;

		// There are no more alignments
		} else { 
			more_alignments = false;
			new_bundle = true; 			// create a fake new bundle start (end of last bundle)
		}

		// There is a new bundle / chromosome
		if (new_bundle || new_chr) {
			hashread.Clear();

			// Process reads in previous bundle
			if (bundle->reads.Count() > 0) {
				bundle->getReady(curr_bundle_start, curr_bundle_end);

				// For -S option: load genome FASTA
				if (use_fasta) {
					const char* chr_name = bundle->refseq.chars();
					GFaSeqGet* fasta_seq = gfasta->fetch(chr_name);
					
					// If failed, try alternative naming convention
					if (fasta_seq == NULL) {
						GStr alt_chr_name;
						
						// If original starts with "chr", try without "chr"
						if (strncmp(chr_name, "chr", 3) == 0) {
							alt_chr_name = chr_name + 3;  // Skip "chr" prefix
						}
						// If original doesn't start with "chr", try adding "chr"
						else {
							alt_chr_name = "chr";
							alt_chr_name.append(chr_name);
						}
						
						GMessage("Trying alternative chromosome name: %s\n", alt_chr_name.chars());
						fasta_seq = gfasta->fetch(alt_chr_name.chars());
					}
					
					if (fasta_seq == NULL) {
						GError("Error: could not retrieve sequence data for %s!\n", 
							bundle->refseq.chars());
					}
					
					bundle->gseq = fasta_seq->copyRange(bundle->start, bundle->end, false, true);
				}
	
				// Process a bunch of bundles in a bundle queue
#ifndef NOTHREADS
				//push this in the bundle queue where it'll be picked up by the threads
				
				int queue_count = 0;

				// Queue Mutex
				queue_mutex.lock();
				bundle_queue->Push(bundle);
				bundle_work |= 0x02; // set bit 1 to 1
				queue_count = bundle_queue->Count();
				queue_mutex.unlock();
				
				// Wait Mutex (wait for a thread to pop this bundle from the queue)
				wait_mutex.lock();
				while (threads_waiting==0) {
					haveThreads.wait(wait_mutex);
				}
				wait_mutex.unlock();
				have_bundles.notify_one();
				
				current_thread::yield();

				// Queue Mutex
				queue_mutex.lock();
				while (bundle_queue->Count()==queue_count) {
					
					queue_mutex.unlock();
					have_bundles.notify_one();
					current_thread::yield();
					queue_mutex.lock();
				}
				queue_mutex.unlock();
				
				// Not using threads, so just process single bundle
#else
				process_bundle(bundle, io);
#endif
			}

			// Clear bundle (no more alignments)
			else { 
				if (THREADING_ENABLED) data_mutex.lock();

				bundle->Clear();

				if (THREADING_ENABLED) {
					clear_data_pool.Push(bundle->idx);
					data_mutex.unlock();
				}
			} 

			// If the chromosome in the BAM file has changed
			if (new_chr) {
			
				// Add guides from chromosome to 'guides'
				total_guides = 0;
				guides = NULL;
				first_possible_overlap = 0;
				// NOTE: What is this below?
				last_guide_idx =- 1;
				
				// Check if BAM ref_id found in guide annotations
				if (refguides.Count() > ref_id) {

					// Check if refguides[ref_id] has guides
					if (refguides[ref_id].rnas.Count() > 0) {
						guides = &(refguides[ref_id].rnas);
						total_guides = guides->Count();
					} else {
						GMessage("Warning: No guides for ref_id=%d (%s), but ref present in guide index.\n",
								ref_id, ref_name);
					}
				} else {
					if (verbose)
						GMessage("Warning: ref_id=%d (%s) not found in guide annotations (refguides.Count() = %d)\n",
								ref_id, ref_name, refguides.Count());
				}
				
				last_ref_name = ref_name;
				last_ref_id = ref_id;
				curr_bundle_end = 0;
			}

			if (!more_alignments) {
				wait_for_bundles();
				break;
			}

			if (THREADING_ENABLED) {
				int new_bidx = wait_for_data(bundles);
				if (new_bidx < 0) {
					//should never happen!
					GError("Error: wait_for_data() returned invalid bundle index(%d)!\n",new_bidx);
					break;
				}
				bundle=&(bundles[new_bidx]);
			}
			
			curr_bundle_start = read_start_pos;
			curr_bundle_end = brec->end;

			// *.** Add guides to bundle
			//GMessage("Total guides: %d\n", total_guides);

			// Move to the first guide that could possibly overlap the read
			first_possible_overlap = last_guide_idx + 1;
			while (first_possible_overlap < total_guides && 
				(int)(*guides)[first_possible_overlap]->end < read_start_pos) {
				first_possible_overlap++; // Skip non-overlapping guides
			}

			int guide_idx = first_possible_overlap;

			// Scan forward to find all guides overlapping the current read
			while (guide_idx < total_guides && 
				(int)(*guides)[guide_idx]->start <= curr_bundle_end) {

				// Expand the current bundle range to include overlapping guides
				curr_bundle_start = std::min(curr_bundle_start, (*guides)[guide_idx]->start);
				curr_bundle_end   = std::max(curr_bundle_end,   (*guides)[guide_idx]->end);

				// On first overlapping guide, check for transitive overlaps *backward*
				if (guide_idx == first_possible_overlap && guide_idx > 0) {
					int back_idx = guide_idx;
					int first_transitive_overlap = guide_idx;

					// Look backward for any guide that also overlaps the updated region
					while (back_idx > last_guide_idx + 1) {
						--back_idx;
						if (curr_bundle_start <= (*guides)[back_idx]->end) {
							first_transitive_overlap = back_idx;
							curr_bundle_start = std::min(curr_bundle_start, (*guides)[back_idx]->start);
						}
					}

					// Add all overlapping and transitively overlapping guides
					for (int i = first_transitive_overlap; i <= guide_idx; ++i) {
						bundle->keepGuide((*guides)[i], &guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
					}
				} else {
					// Directly overlapping guide
					bundle->keepGuide((*guides)[guide_idx], &guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
				}
				++guide_idx;
			}

			// Update last processed guide index
			last_guide_idx = guide_idx - 1;

			bundle->refseq = last_ref_name;
			bundle->start = curr_bundle_start;
			bundle->end = curr_bundle_end;
		} 

		

		// Current read extends the bundle
		// This might not happen if a longer guide had already been added to the bundle
		if (curr_bundle_end < (int)brec->end) {
			curr_bundle_end = brec->end;

			// Add newly overlapping guides to the bundle
			// Keep expanding bundle to the right as long as new guides overlap with the expanded region
			while (true) {
				bool expanded = false;
				while (last_guide_idx + 1 < total_guides && (int)(*guides)[last_guide_idx + 1]->start <= curr_bundle_end) {
					++last_guide_idx;
					auto* guide = (*guides)[last_guide_idx];
					if (guide->end >= curr_bundle_start) {
						bundle->keepGuide(guide, &guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
						if (curr_bundle_end < guide->end) {
							curr_bundle_end = guide->end;
							expanded = true;
						}
					}
				}
				if (!expanded) break;
			}
		} 

		GReadAlnData alndata(brec, 0, nh, hi);
     	bool overlaps_guide = bundle->evalReadAln(alndata, splice_strand);
		
		// Reads with "." strand never overlap
		// Only process the read if it overlaps a guide
		if (overlaps_guide) {

			// Splice_strand may have been set by evalReadAln
			if (splice_strand == '+') 
				alndata.strand = 1;
			else if (splice_strand == '-') 
				alndata.strand = -1;

			process_read_in(curr_bundle_start, curr_bundle_end, *bundle, hashread, alndata, brec);
		}
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Finish and clean up
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// Delete thread and bundle arrays
	if (THREADING_ENABLED) {
		for (int t = 0; t < num_cpus; t++) threads[t].join();
		if (verbose) {
			printTime(stderr);
			GMessage(" All threads finished.\n");
		}
		delete[] threads;
		delete[] bundles;
		delete bundle_queue;

	}

	io->stop(); 		// close all BAM files
	delete io;
	delete worker_args;
	delete gfasta;
	
	// Print time, declare finished
	if (verbose) {
		printTime(stderr);
		GMessage(" Done.\n");
	}

	//if (c_out && c_out!=stdout) fclose(c_out);  
	// eite

	f_out=stdout;
	// if(outfname!="stdout") {
	// 	f_out=fopen(outfname.chars(), "w");
	// 	if (f_out==NULL) GError("Error creating output file %s\n", outfname.chars());
	// }

	fprintf(f_out,"# ");
	// args.printCmdLine(f_out);
	fprintf(f_out,"# Bramble version %s\n",VERSION);

	//FILE *g_out=NULL;
	//FILE* tmp_fin=fopen(tmp_fname.chars(),"rt");
	//fclose(f_out);
	gffnames_unref(guide_seq_names); //deallocate names collection

	if (GMEMTRACE && verbose) GMessage(" Max bundle memory: %6.1fMB for bundle %s\n", maxMemRS/1024, maxMemBundle.chars());

	// DONT close f again it's already closed somewhere
}
