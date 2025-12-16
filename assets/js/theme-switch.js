window.addEventListener("DOMContentLoaded", function() {
    const toggleDarkMode = document.getElementById("theme-toggle");
  const themeStyleLink = document.getElementById('theme-style');

    // If the toggle isn't present on the page, nothing to do
    if (!toggleDarkMode) return;

    if (localStorage.getItem('theme') === 'dark') {
      setTheme('dark');
    } else {
      setTheme('light');
    }

    // attach click handler using jtd helper if available, otherwise fallback
    if (typeof jtd !== 'undefined' && typeof jtd.addEvent === 'function') {
      jtd.addEvent(toggleDarkMode, 'click', function(){
        const currentTheme = getTheme();
        const newTheme = currentTheme === 'dark' ? 'light' : 'dark';

        localStorage.setItem('theme', newTheme);
        setTheme(newTheme);
        if (typeof jtd.setTheme === 'function') jtd.setTheme(newTheme);
      });
    } else {
      toggleDarkMode.addEventListener('click', function(){
        const currentTheme = getTheme();
        const newTheme = currentTheme === 'dark' ? 'light' : 'dark';

        localStorage.setItem('theme', newTheme);
        setTheme(newTheme);
      });
    }
  
    function getTheme() {
      return document.documentElement.classList.contains('dark-mode') ? 'dark' : 'light';
    }
  
    function setTheme(theme) {
      if (theme === 'dark') {
        toggleDarkMode.innerHTML = `<svg width='18px' height='18px'><use href="#svg-moon"></use></svg>`;
        document.documentElement.classList.add('dark-mode');
        document.documentElement.classList.remove('light-mode');
        if (themeStyleLink) {
          // Prefer dataset values injected by Liquid in the head
          if (themeStyleLink.dataset && themeStyleLink.dataset.dark) {
            themeStyleLink.href = themeStyleLink.dataset.dark;
          } else {
            themeStyleLink.href = '/assets/css/just-the-docs-dark.css';
          }
        }
      } else {
        toggleDarkMode.innerHTML = `<svg width='18px' height='18px'><use href="#svg-sun"></use></svg>`;
        document.documentElement.classList.add('light-mode');
        document.documentElement.classList.remove('dark-mode');
        if (themeStyleLink) {
          if (themeStyleLink.dataset && themeStyleLink.dataset.light) {
            themeStyleLink.href = themeStyleLink.dataset.light;
          } else {
            themeStyleLink.href = '/assets/css/just-the-docs-light.css';
          }
        }
      }
    }
  });