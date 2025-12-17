document.addEventListener("DOMContentLoaded", function () {
  if (document.querySelector('#toc')) {
    tocbot.init({
      tocSelector: '#toc',
      contentSelector: '.main-content',
      headingSelector: 'h2, h3, h4',
      collapseDepth: 2,
      scrollSmooth: true
    });
  }
});
