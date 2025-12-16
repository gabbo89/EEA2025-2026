window.addEventListener('DOMContentLoaded', function() {
  // find all potential theme toggle elements (sidebar or header)
  const toggles = Array.from(document.querySelectorAll('#theme-toggle, #theme-toggle-header, .theme-toggle'));
  const themeStyleLink = document.getElementById('theme-style');

  if (!toggles.length) return; // nothing to do

  const sunSVG = `<svg width='18px' height='18px'><use href="#svg-sun"></use></svg>`;
  const moonSVG = `<svg width='18px' height='18px'><use href="#svg-moon"></use></svg>`;

  function getTheme() {
    return document.documentElement.classList.contains('dark-mode') ? 'dark' : 'light';
  }

  function setTheme(theme) {
    // update all toggle icons
    toggles.forEach(el => {
      try { el.innerHTML = theme === 'dark' ? moonSVG : sunSVG; } catch (e) { /* ignore DOM errors */ }
    });

    if (theme === 'dark') {
      document.documentElement.classList.add('dark-mode');
      document.documentElement.classList.remove('light-mode');
      if (themeStyleLink && themeStyleLink.dataset && themeStyleLink.dataset.dark) {
        themeStyleLink.href = themeStyleLink.dataset.dark;
      }
    } else {
      document.documentElement.classList.add('light-mode');
      document.documentElement.classList.remove('dark-mode');
      if (themeStyleLink && themeStyleLink.dataset && themeStyleLink.dataset.light) {
        themeStyleLink.href = themeStyleLink.dataset.light;
      }
    }
  }

  // initialize theme from localStorage or system preference
  const initial = localStorage.getItem('theme') || (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
  if (!localStorage.getItem('theme')) localStorage.setItem('theme', initial);
  setTheme(initial);

  // attach click handlers
  toggles.forEach(el => {
    const handler = function() {
      const newTheme = getTheme() === 'dark' ? 'light' : 'dark';
      localStorage.setItem('theme', newTheme);
      setTheme(newTheme);
      if (typeof jtd !== 'undefined' && typeof jtd.setTheme === 'function') jtd.setTheme(newTheme);
    };

    if (typeof jtd !== 'undefined' && typeof jtd.addEvent === 'function') {
      try { jtd.addEvent(el, 'click', handler); } catch (e) { el.addEventListener('click', handler); }
    } else {
      el.addEventListener('click', handler);
    }
  });
});