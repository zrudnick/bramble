// Site navigation: the one place to add, rename, or reorder pages.
// Each page has <nav id="sidebar"></nav> and <div id="pager"></div>; this file
// fills them in. The current page's section starts open, the others closed.
const NAV = [
  {
    section: "Getting started",
    pages: [
      ["index.html", "Introduction"],
      ["installation.html", "Installation"],
      ["quick-start.html", "Quick start"],
    ],
  },
  {
    section: "Command-line interface",
    pages: [
      ["cli-required.html", "Required arguments"],
      ["cli-presets.html", "Read-type presets"],
      ["cli-boundaries.html", "Boundary tolerances"],
      ["cli-long-reads.html", "Long-read options"],
      ["cli-clip-rescue.html", "Clip rescue"],
      ["cli-strandedness.html", "Strandedness"],
      ["cli-other.html", "Threads and other options"],
    ],
  },
  {
    section: "Output",
    pages: [["output-bam.html", "Transcriptome BAM"]],
  },
  {
    section: "Guides",
    pages: [["guides.html", "About guides"]],
  },
];

(function () {
  const here = location.pathname.split("/").pop() || "index.html";
  const flat = NAV.flatMap((s) => s.pages);

  const sidebar = document.getElementById("sidebar");
  if (sidebar) {
    for (const s of NAV) {
      const details = document.createElement("details");
      details.open = s.pages.some(([href]) => href === here);
      const summary = document.createElement("summary");
      summary.textContent = s.section;
      details.appendChild(summary);
      const ul = document.createElement("ul");
      for (const [href, title] of s.pages) {
        const li = document.createElement("li");
        const a = document.createElement("a");
        a.href = href;
        a.textContent = title;
        if (href === here) a.className = "active";
        li.appendChild(a);
        ul.appendChild(li);
      }
      details.appendChild(ul);
      sidebar.appendChild(details);
    }
  }

  // previous / next links, in sidebar order
  const pager = document.getElementById("pager");
  const i = flat.findIndex(([href]) => href === here);
  if (pager && i >= 0) {
    if (i > 0) {
      const a = document.createElement("a");
      a.href = flat[i - 1][0];
      a.className = "prev";
      a.textContent = "← " + flat[i - 1][1];
      pager.appendChild(a);
    }
    if (i < flat.length - 1) {
      const a = document.createElement("a");
      a.href = flat[i + 1][0];
      a.className = "next";
      a.textContent = flat[i + 1][1] + " →";
      pager.appendChild(a);
    }
  }

  // small screens: the menu button shows and hides the sidebar
  const button = document.getElementById("menu-button");
  if (button && sidebar) {
    button.addEventListener("click", () => {
      const open = sidebar.classList.toggle("open");
      button.setAttribute("aria-expanded", open);
    });
  }
})();
