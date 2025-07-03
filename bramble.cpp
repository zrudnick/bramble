// Bramble v1.0.0
// Main file

#include <numeric>

#include "bramble.h"
#include "bam.h"
#include "reads.h"
#include "tree.h"
#include "GThreads.h" 			// for THREADING_ENABLED
#include "htslib/sam.h"

#define VERSION "1.0.0"

#define USAGE "Bramble v" VERSION " usage:\n\n\
bramble <in.bam ..> [-G <guide_gff>] [-o <out.bam>] [-p <cpus>] [-S <genome.fa>] \n\
 [--help] [--version] [--verbose] [--long] [--fr] [--rf] \n\
 \n\
Convert BAM file coordinates from genomic to transcriptomic.\n\
Options:\n\
 --help     : print this usage message and exit\n\
 --version  : print just the version at stdout and exit\n\
 --verbose  : verbose (log bundle processing details)\n\
 --long     : BAM file contains long reads\n\
 --fr       : assume stranded library fw-firststrand\n\
 --rf       : assume stranded library fw-secondstrand\n\
 -G <file>  : reference annotation to use for guiding the BAM conversion (GTF/GFF)\n\
 -o <file>  : output path/file name for the converted BAM file (default: stdout)\n\
 -p <int>   : number of threads (CPUs) to use (default: 1)\n\
 -S <file>  : genome sequence file (FASTA format)\n\
"

FILE* F_OUT=NULL;
bool VERBOSE = false;				// Verbose, --verbose
bool LONG_READS = false;				// BAM file contains long reads, --long
bool FR_STRAND = false;				// Read 1 is on forward strand, --fr
bool RF_STRAND = false;				// Read 1 is on reverse strand, --fr
uint8_t N_THREADS = 1;				// Threads, -p
GStr BAM_PATH_OUT;					// Output BAM path, -o
GStr GUIDE_GFF; 					// Reference annotation, -G
GFastaDb* GFASTA = NULL;			// FASTA file, -S
bool USE_FASTA = false;				// Use FASTA for reducing soft clips
const char* BAM_FILE_IN;			// Input BAM path
const char* BAM_FILE_OUT;			// Output BAM path
const char* OUT_DIR;				// Output folder

#ifndef NOTHREADS
// Threading: single producer, multiple consumers
// Main thread/program is always loading the producer

// THREADING_ENABLED defined in bramble.h

GMutex DATA_MUTEX; 					// Manage availability of data records ready to be loaded by main thread
GVec<int> CLEAR_DATA_POOL; 			// Indices of data bundles cleared for loading by main thread (clear data pool)
GConditionVar HAVE_BUNDLES; 		// Will notify a thread that a bundle was loaded in the ready queue
                           			// (or that no more bundles are coming)
char BUNDLE_WORK = 1; 				// Bit 0 set if bundles are still being prepared (BAM file not exhausted yet)
                  					// Bit 1 set if there are Bundles ready in the queue

GMutex WAIT_MUTEX; 					// Controls THREADS_WAITING (idle threads counter)
uint8_t THREADS_WAITING; 			// Idle worker threads
GConditionVar HAVE_THREADS; 		// Will notify the bundle loader when a thread
                          			// Is available to process the currently loaded bundle

GConditionVar HAVE_CLEAR; 			// Will notify when bundle buf space available
GMutex QUEUE_MUTEX; 				// Controls bundle_queue and bundles access
GFastMutex LOG_MUTEX; 				// Only when verbose - to avoid mangling the log output
GFastMutex READING_MUTEX;
#endif

bool NO_MORE_BUNDLES = false;

// Process input parameters
void process_options(GArgs& args) {
	// Help
	if (args.getOpt("help")) {
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
	if (args.getOpt("fr")) {
		FR_STRAND = true;
	}

	// RF (Reverse-Forward)
	// First read is on reverse strand, second read is on forward strand
	if (args.getOpt("rf")) {
		RF_STRAND = true;
		if (FR_STRAND) {
			GError("Error: --fr and --rf options are incompatible.\n");
		}
	}

	// Verbose
	if (args.getOpt("verbose")) {
		VERBOSE = true;
		fprintf(stderr, "Running Bramble " VERSION ". Command line:\n");
		args.printCmdLine(stderr);
	}

	// Long Reads
    if (args.getOpt("long")) {
        LONG_READS = true;
    }

	// Number of Threads, -p
	GStr s = args.getOpt('p');
	if (!s.is_empty()) {
		N_THREADS = s.asInt();
		if (N_THREADS <= 0) N_THREADS = 1;
	}

	// Reference Annotation, -G
    s = args.getOpt('G');
    if (!s.is_empty()) {
        GUIDE_GFF = s;
        if (!fileExists(GUIDE_GFF.chars()) > 1) {
            GError("Error: reference annotation file (%s) not found.\n",
                   GUIDE_GFF.chars());
        }
    }

	// Genome sequence file, -S
	s = args.getOpt('S');
	if (!s.is_empty()) {
		GFASTA = new GFastaDb(s.chars());
		USE_FASTA = true;
	}

	// Output path, -o
    s = args.getOpt('o');
    if (!s.is_empty()) {
        BAM_PATH_OUT = s;
		BAM_FILE_OUT = BAM_PATH_OUT;
        GMessage("Output file: %s\n", BAM_FILE_OUT);
    } else {
        GMessage("%s\nError: no output file name provided!\n", USAGE);
        exit(1);
    }


	// Validate input files
	int n_bam = args.startNonOpt();
	if (n_bam < 1) {
		GMessage("%s\nError: no input file provided!\n", USAGE);
		exit(1);
	}

	if (GUIDE_GFF == NULL) {
		GMessage("%s\nError: no input file provided!\n", USAGE);
		exit(1);
	}

	// Process input alignment files
	const char* ifn = NULL;
	int i = 0;
	while ((ifn = args.nextNonOpt()) != NULL) {
		if (i == 0) BAM_FILE_IN = ifn;				// currently only uses first BAM file
		i++;
	}
}

// Check if there are more bundles
bool has_more_bundles() {
	if (THREADING_ENABLED) GLockGuard<GFastMutex> lock(READING_MUTEX);
	return !NO_MORE_BUNDLES;
}

// Wait for threads and set bundle status
void wait_for_bundles() {

	if (THREADING_ENABLED) {
		READING_MUTEX.lock();
		NO_MORE_BUNDLES = true;
		READING_MUTEX.unlock();

		QUEUE_MUTEX.lock();
		BUNDLE_WORK &= ~(int)0x01; // clear bit 0;
		QUEUE_MUTEX.unlock();

		bool are_threads_waiting = true;
		while (are_threads_waiting) {

			WAIT_MUTEX.lock();
			are_threads_waiting = (THREADS_WAITING > 0);
			WAIT_MUTEX.unlock();

			if (are_threads_waiting) {
				HAVE_BUNDLES.notify_all();
				current_thread::sleep_for(1);

				WAIT_MUTEX.lock();
				are_threads_waiting = (THREADS_WAITING > 0);
				WAIT_MUTEX.unlock();

				current_thread::sleep_for(1);
			}
		}
	} 
	
	else {
		NO_MORE_BUNDLES = true;
	}
}

// Process current bundle
void process_bundle(BundleData* bundle, BamIO* io) {

	if (VERBOSE) {
		if (THREADING_ENABLED) GLockGuard<GFastMutex> lock(LOG_MUTEX);
		GMessage(">bundle %s:%d-%d [%lu alignments (%d distinct), %d junctions, %d guides] begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->reads.Count(), bundle->junction.Count(),
				bundle->guides.Count());
	}

	// This is where reads are added to new BAM file
	convert_reads(bundle, io);

	if (VERBOSE) {
		if (THREADING_ENABLED) GLockGuard<GFastMutex> lock(LOG_MUTEX);
		GMessage("^bundle %s:%d-%d done (%d processed potential transcripts).\n", bundle->refseq.chars(), bundle->start, bundle->end);
	}

	bundle->Clear();
}

// Check that there aren't any threads waiting
bool no_threads_waiting() {
	WAIT_MUTEX.lock();
	int threads = THREADS_WAITING;
	WAIT_MUTEX.unlock();
	return (threads < 1);
}

// Worker thread waits for incoming bundle, calls process_bundle
void worker_thread(GThreadData& td) {
	WorkerArgs* args = (WorkerArgs*)(td.udata);
	GPVec<BundleData>* bundle_queue = args->bundle_queue;
	BamIO* io = args->io;

	// Wait for a ready bundle in the queue
	QUEUE_MUTEX.lock(); // enter wait-for-notification loop

	while (BUNDLE_WORK) {
		WAIT_MUTEX.lock();
		THREADS_WAITING++;
		QUEUE_MUTEX.unlock();
		WAIT_MUTEX.unlock();
		HAVE_THREADS.notify_one(); 			// in case main thread is waiting
		current_thread::yield();
		QUEUE_MUTEX.lock();
		while (BUNDLE_WORK && bundle_queue->Count() == 0) {
			// Unlocks QUEUE_MUTEX and wait until notified
			// When notified, locks QUEUE_MUTEX and resume
			HAVE_BUNDLES.wait(QUEUE_MUTEX); 
		}

		WAIT_MUTEX.lock();
		if (THREADS_WAITING > 0) THREADS_WAITING--;
		WAIT_MUTEX.unlock();

		BundleData* readyBundle = NULL;
		if ((BUNDLE_WORK & 0x02) != 0) {
			readyBundle = bundle_queue->Pop();

			if (readyBundle != NULL) {
				if (bundle_queue->Count() == 0)
				BUNDLE_WORK &= ~(int)0x02; 		// clear bit 1 (queue is empty)
			
				QUEUE_MUTEX.unlock();
				process_bundle(readyBundle, io);

				DATA_MUTEX.lock();
				CLEAR_DATA_POOL.Push(readyBundle->idx);
				DATA_MUTEX.unlock();
				
				HAVE_CLEAR.notify_one(); 		// inform main thread
				current_thread::yield();
				
				QUEUE_MUTEX.lock();
			}
		}
	}
	QUEUE_MUTEX.unlock();
}

// Prepare the next available bundle slot for loading
int wait_for_data(BundleData* bundles) {
	int idx = -1;

	DATA_MUTEX.lock();
	while (CLEAR_DATA_POOL.Count() == 0) {
		HAVE_CLEAR.wait(DATA_MUTEX);
	}
	idx = CLEAR_DATA_POOL.Pop();
	if (idx >= 0) bundles[idx].status = BundleStatus::BUNDLE_STATUS_LOADING;
	DATA_MUTEX.unlock();

	return idx;
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

int main(int argc, char* argv[]) {

	GArgs args(argc, argv,
        "help;version;fr;rf;verbose;long;p:G:S:o:");
 	args.printError(USAGE, true);
 	process_options(args);

	const char* ERR_BAM_SORT = "\nError: the input alignment file is not sorted!\n";

	// Table indexes for raw counts data
	GPVec<RC_TData> guides_RC_tdata(true);     	 		// Raw count data for all guide transcripts
	GPVec<RC_Feature> guides_RC_exons(true);    		// Raw count data for all guide exons
	GPVec<RC_Feature> guides_RC_introns(true);  		// Raw count data for all guide introns

	// ^**^^^*^**^^^*^**^^^*^**^^^*
	// Read in reference annotation
	// ^**^^^*^**^^^*^**^^^*^**^^^*

	// Create SAM header file
	const GStr header_path = "tx_header.sam";
    FILE* header_file = fopen(header_path.chars(), "w");
    if (header_file == NULL) GError("Error creating file: %s\n", header_path.chars());
    fprintf(header_file, "@HD\tVN:1.0\tSO:coordinate\n");
	
	if (VERBOSE) GMessage(" Loading reference annotation (guides)..\n");
	
	// Open GFF/GTF file
	FILE* f = fopen(GUIDE_GFF.chars(), "r");
	if (f == NULL) GError("Error: could not open reference annotation file (%s)!\n", GUIDE_GFF.chars());

	// GffReader: transcripts only, sort by location
	GffReader gffreader(f, true, true);           // loading only recognizable transcript features
	gffreader.setRefAlphaSorted();                // alphabetical sorting of RefSeq IDs
	gffreader.showWarnings(VERBOSE);

	// keep attributes, merge close exons, no exon attributes
	// merge_close_exons must be false for correct transcriptome header construction!
	gffreader.readAll(true, false, false);
	
	int n_refguides = gffreader.gseqtable.Count();
	if (n_refguides == 0 || gffreader.gflst.Count() == 0) {
		GError("Error: could not any valid reference transcripts in %s (invalid GTF/GFF file?)\n", GUIDE_GFF.chars());
	}

	GVec<GRefData> refguides;           // vector with guides for each chromosome
	refguides.setCount(n_refguides); 	// maximum reference guide ID

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
		if (guide->exons.Count() == 0) {
			if (VERBOSE)
				GMessage("Warning: exonless GFF %s feature with ID %s found, added implicit exon %d-%d.\n",
							guide->getFeatureName(), guide->getID(), guide->start, guide->end);
			guide->addExon(guide->start, guide->end); // should never happen!
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
	   	grefdata.add(&gffreader, guide); // transcripts already sorted by location

	}

	GffNames* guide_seq_names = GffObj::names; 		// might have been populated already by gff data
	gffnames_ref(guide_seq_names);  		// initialize the names collection if not guided

	fprintf(header_file, "@CO\tGenerated from GTF: %s\n", GUIDE_GFF.chars());
    fclose(header_file);

	if (VERBOSE) {
		GMessage("Transcriptome header created\n");
		GMessage(" %d reference transcripts loaded\n", gffreader.gflst.Count());
	}

	// ^**^^^*^**^^^*^**^^^*^**^^^*
	// Start BAM IO
	// ^**^^^*^**^^^*^**^^^*^**^^^*

	BamIO* io = new BamIO(BAM_FILE_IN, BAM_FILE_OUT, header_path);
	io->start();

	// *.** Bundle reads into overlapping groups

	GHash<int> hashread;      			// read_name:pos:hit_index => readlist index
	GList<GffObj>* guides = NULL; 		// list of transcripts on a specific reference
	uint curr_bundle_start = 0;
	uint curr_bundle_end = 0;
	int first_possible_overlap = 0;
	int last_guide_idx = -1;
	int total_guides = 0;

	GStr last_ref_name;
	int last_ref_id = -1; 			// last seen ref_id

#ifndef NOTHREADS // THREADING_ENABLED
	#define DEF_TSTACK_SIZE 8388608
	size_t def_stack_size = DEF_TSTACK_SIZE;

	#ifdef _GTHREADS_POSIX_
	size_t tstack_size = GThread::defaultStackSize();
	if (tstack_size < DEF_TSTACK_SIZE) def_stack_size = DEF_TSTACK_SIZE;
	if (VERBOSE) {
		if (tstack_size < def_stack_size) {
			GMessage("Default stack size for threads: %d (increased to %d)\n", 
				tstack_size, def_stack_size);
		}
		else GMessage("Default stack size for threads: %d\n", tstack_size);
	}
	#endif

	GThread* threads = new GThread[N_THREADS]; 							// threads for processing bundles
	GPVec<BundleData>* bundle_queue = new GPVec<BundleData>(false); 	// queue for holding loaded bundles
	BundleData* bundles = new BundleData[N_THREADS + 1]; 				// redef with more bundles

	CLEAR_DATA_POOL.setCapacity(N_THREADS + 1);

	WorkerArgs* worker_args = new WorkerArgs(bundle_queue, io);

	// Start worker threads
	for (int b = 0; b < N_THREADS; b++) {
		threads[b].kickStart(worker_thread, (void*) worker_args, def_stack_size);
		bundles[b+1].idx = b + 1;
		CLEAR_DATA_POOL.Push(b);
	}
	BundleData* bundle = &(bundles[N_THREADS]);

	
#else // !THREADING_ENABLED

	// Just put everything into the same bundle
	BundleData bundles[1];
	BundleData* bundle = &(bundles[0]);
#endif
	
	GSamRecord* brec = nullptr;
	bool more_alignments = true;
	int prev_pos = 0;

	while (more_alignments) {
		const char* ref_name = NULL;
		char splice_strand = 0;
		int read_start_pos = 0;
		int nh = 1;
		int hi = 0;
		int ref_id = last_ref_id;  // current chromosome ID

		bool new_bundle = false;
		bool new_chromosome = false;
		
		if ((brec = io->next()) != NULL) {
			if (brec->isUnmapped()) continue;

			ref_name = brec->refName();
			read_start_pos = brec->start;
			splice_strand = brec->spliceStrand(); 	// tagged strand gets priority
			
			// Set strand if stranded library
			if ((splice_strand == '.') && (FR_STRAND || RF_STRAND)) { 

				if (brec->isPaired()) {
					if (brec->pairOrder() == 1) { 	// first read in pair
						if ((RF_STRAND && brec->revStrand()) || (FR_STRAND && !brec->revStrand())) 
							 splice_strand = '+';
						else splice_strand = '-';
					}
					else {					    	// second read in pair
						if ((RF_STRAND && brec->revStrand())||(FR_STRAND && !brec->revStrand())) 
							 splice_strand = '-';
						else splice_strand = '+';
					}

				// Read isn't paired
				} else {
					if ((RF_STRAND && brec->revStrand()) || (FR_STRAND && !brec->revStrand())) 
						 splice_strand = '+';
					else splice_strand = '-';
				}

			}

			if (ref_name == NULL) {
				GError("Error: cannot retrieve target seq name from BAM record!\n");
			}

			// Are we at a new chromosome?
			new_chromosome = (last_ref_name.is_empty() || last_ref_name != ref_name);
			if (new_chromosome) {
				ref_id = guide_seq_names->gseqs.addName(ref_name);
				
				if (ref_id >= n_refguides) {
					if (VERBOSE) {
						GMessage("WARNING: no reference transcripts found for genomic sequence \"%s\"!\n",
								ref_name); 
						}
				}
				
				prev_pos = 0;
			}
			if (read_start_pos < prev_pos) GError("%s\nread %s (start %d) found at position %d on %s when prev_pos=%d\n",
				ERR_BAM_SORT, brec->name(), brec->start, read_start_pos, ref_name, prev_pos);
			prev_pos = read_start_pos;

			nh = brec->tag_int("NH");
			nh =  nh ? nh != 0 : 1;  		// number of hits
			hi = brec->tag_int("HI");		// query hit index

			if ((!new_chromosome) && (curr_bundle_end > 0) && 
			    (read_start_pos > (curr_bundle_end + (int)RUNOFF_DIST))) {
					new_bundle = true;
			}

		// No more alignments
		} else { 
			more_alignments = false;
			new_bundle = true; 			// if last bundle, create fake start
		}

		// *.** New bundle / chromosome

		if (new_bundle || new_chromosome) {
			hashread.Clear();

			// Process reads in previous bundle
			if (bundle->reads.Count() > 0) {
				bundle->getReady(curr_bundle_start, curr_bundle_end);

				// For -S option: load genome FASTA
				if (USE_FASTA) {
					const char* chr_name = bundle->refseq.chars();
					GFaSeqGet* fasta_seq = GFASTA->fetch(chr_name);
					
					// Try alternative naming convention
					if (fasta_seq == NULL) {
						GStr alt_chr_name;
						if (strncmp(chr_name, "chr", 3) == 0) {
							alt_chr_name = chr_name + 3;  // skip "chr" prefix
						}
						else {
							alt_chr_name = "chr";
							alt_chr_name.append(chr_name); // add "chr" prefix
						}
						fasta_seq = GFASTA->fetch(alt_chr_name.chars());
					}
					
					if (fasta_seq == NULL) {
						GError("Error: could not retrieve sequence data for %s!\n", 
							bundle->refseq.chars());
					}
					
					bundle->gseq = fasta_seq->copyRange(bundle->start, bundle->end, false, true);
				}
	
#ifndef NOTHREADS // THREADING_ENABLED
				
				// Push this in the bundle queue where it'll be picked up by the threads	
				int queue_count = 0;

				QUEUE_MUTEX.lock();
				bundle_queue->Push(bundle);
				BUNDLE_WORK |= 0x02; // set bit 1 to 1
				queue_count = bundle_queue->Count();
				QUEUE_MUTEX.unlock();
				
				WAIT_MUTEX.lock();
				while (THREADS_WAITING == 0) {
					HAVE_THREADS.wait(WAIT_MUTEX);
				}
				WAIT_MUTEX.unlock();
				HAVE_BUNDLES.notify_one();
				
				current_thread::yield();

				QUEUE_MUTEX.lock();
				while (bundle_queue->Count()==queue_count) {
					QUEUE_MUTEX.unlock();
					HAVE_BUNDLES.notify_one();
					current_thread::yield();
					QUEUE_MUTEX.lock();
				}
				QUEUE_MUTEX.unlock();
				
				
#else // !THREADING_ENABLED
				
				// Just process single bundle
				process_bundle(bundle, io);
#endif
			}

			// Clear bundle (no more alignments)
			else { 
				if (THREADING_ENABLED) DATA_MUTEX.lock();

				bundle->Clear();

				if (THREADING_ENABLED) {
					CLEAR_DATA_POOL.Push(bundle->idx);
					DATA_MUTEX.unlock();
				}
			} 

			// If there is a new chromosome
			if (new_chromosome) {
			
				// Add guides from chromosome
				total_guides = 0;
				guides = NULL;
				first_possible_overlap = 0;
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
					if (VERBOSE) {
						GMessage("Warning: ref_id=%d (%s) not found in guide annotations (refguides.Count() = %d)\n",
								ref_id, ref_name, refguides.Count());
					}	
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
					GError("Error: wait_for_data() returned invalid bundle index(%d)!\n",new_bidx);
					break; // should never happen!
				}
				bundle = &(bundles[new_bidx]);
			}
			
			curr_bundle_start = read_start_pos;
			curr_bundle_end = brec->end;

			// *.** Add guides to bundle

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
		
		// Only process the read if it overlaps a guide
		if (overlaps_guide) { // reads with "." strand never overlap

			// Splice_strand may have been set by evalReadAln
			if (splice_strand == '+') 
				alndata.strand = 1;
			else if (splice_strand == '-') 
				alndata.strand = -1;

			process_read_in(curr_bundle_start, curr_bundle_end, *bundle, hashread, alndata, brec);
		}
	}

	//^**^^^*^**^^^*^**^^^*^**^^^*
	// Finish and clean up
	//^**^^^*^**^^^*^**^^^*^**^^^*

	// Delete thread and bundle arrays
	if (THREADING_ENABLED) {
		for (int t = 0; t < N_THREADS; t++) threads[t].join();
		if (VERBOSE) {
			GMessage(" All threads finished.\n");
		}
		delete[] threads;
		delete[] bundles;
		delete bundle_queue;
	}

	io->stop(); 		// close BAM reader & writer
	delete io;
	delete worker_args;
	delete GFASTA;

	F_OUT = stdout;
	fprintf(F_OUT, "# ");
	args.printCmdLine(F_OUT);
	fprintf(F_OUT,"# Bramble version %s\n", VERSION);

	gffnames_unref(guide_seq_names); 	// deallocate names collection
}