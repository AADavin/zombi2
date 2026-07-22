// MathJax for the display formulas in the chapters (the P(·) factorisations of Chapter 2, the rate
// grammar). pymdownx.arithmatex runs in `generic: true` mode, which hands MathJax \(…\) and \[…\]
// wrapped in .arithmatex, so only those spans are typeset and the rest of the page is left alone.
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex",
  },
};

// Material swaps pages in without a reload, so re-typeset on every navigation.
document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
