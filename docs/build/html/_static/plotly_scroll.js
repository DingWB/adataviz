// Safer Plotly wheel forwarding.
// Attach wheel handlers directly to each .plotly-graph-div (existing and future)
// to avoid global capture side-effects. Forward the event to the inner SVG,
// mark forwarded events to avoid recursion, and briefly focus+blur the SVG so
// Plotly picks up scroll-zoom without permanently trapping input.

(function () {
  const FORWARD_FLAG = '__plotly_forwarded';

  function handleWheelOnPlot(e) {
    // Ignore events that were already forwarded
    if (e[FORWARD_FLAG]) return;

    try { console.debug && console.debug('[plotly_scroll] wheel received', { deltaX: e.deltaX, deltaY: e.deltaY, target: e.target }); } catch (e) {}

    const plot = e.currentTarget;
    // Heuristic: detect UMAP plots by searching nearby text for the word 'UMAP'.
    // For UMAP plots we avoid stealing the wheel unless the user holds Ctrl/Meta
    // (so page scroll remains usable and wheel-zoom requires an explicit modifier).
    function isUMAPPlot(plotEl) {
      try {
        const container = plotEl.closest && (plotEl.closest('.output') || plotEl.closest('.output_area')) || plotEl.parentElement;
        if (!container) return false;
        const txt = (container.innerText || '').slice(0, 400).toLowerCase();
        return txt.indexOf('umap') !== -1;
      } catch (err) { return false; }
    }

    const isUMAP = isUMAPPlot(plot);

    // For UMAP: only forward wheel events when Ctrl/Meta is held. Otherwise
    // allow the wheel to scroll the page (do not preventDefault).
    if (isUMAP && !(e.ctrlKey || e.metaKey)) {
      return;
    }

    // Only act for direct wheel events on the plot container
    e.preventDefault();

    const svg = plot.querySelector('svg');
    if (!svg) return;

    // Build a synthetic wheel event. Allow bubbling so Plotly internal listeners
    // that rely on event propagation receive it. We still mark the event so
    // our handler can ignore it if it's seen again.
    const evt = new WheelEvent('wheel', {
      deltaX: e.deltaX,
      deltaY: e.deltaY,
      deltaMode: e.deltaMode || 0,
      clientX: e.clientX,
      clientY: e.clientY,
      bubbles: true,
      cancelable: true,
      composed: true
    });

    try {
      Object.defineProperty(evt, FORWARD_FLAG, { value: true });
    } catch (err) {
      evt[FORWARD_FLAG] = true;
    }

    // Focus the SVG briefly so Plotly's handlers that check focus will work,
    // then blur after a short delay to avoid trapping keyboard/scroll focus.
    try { svg.focus && svg.focus(); } catch (err) {}

    svg.dispatchEvent(evt);

    // Restore focus to body shortly after to avoid trapping input on the plot.
    setTimeout(() => { try { document.activeElement && document.activeElement.blur && document.activeElement.blur(); } catch (err) {} }, 50);
  }

  // Global capture-phase wheel handler to prevent Plotly from stealing
  // the wheel when the target is inside a UMAP plot and the user is not
  // explicitly holding Ctrl/Meta. This runs in capture and stops immediate
  // propagation so downstream Plotly handlers don't preventDefault.
  function captureWheelForUMAP(e) {
    try {
      if (e.ctrlKey || e.metaKey) return; // user explicitly wants zoom
      var el = e.target;
      while (el && el !== document.documentElement) {
        if (el.classList && el.classList.contains && el.classList.contains('plotly-graph-div')) {
          try {
            // Heuristic: check for nearby 'UMAP' text inside the plot's container
            const container = el.closest && (el.closest('.output') || el.closest('.output_area')) || el.parentElement;
            const txt = (container && (container.innerText || '') || '').toLowerCase();
            if (txt.indexOf('umap') !== -1) {
              // Prevent other listeners (including Plotly) from running so
              // they cannot call preventDefault and trap the page scroll.
              e.stopImmediatePropagation();
              // Do not call preventDefault here — allow the document to scroll.
              return;
            }
          } catch (err) { return; }
        }
        el = el.parentElement;
      }
    } catch (err) {}
  }

  try {
    document.addEventListener('wheel', captureWheelForUMAP, { passive: false, capture: true });
  } catch (err) {}

  function attachToPlot(plot) {
    if (!plot.__plotly_scroll_attached) {
      plot.addEventListener('wheel', handleWheelOnPlot, { passive: false });
      plot.__plotly_scroll_attached = true;
      try { console.debug && console.debug('[plotly_scroll] attached wheel handler', plot); } catch (e) {}
    }
  }

  // Ensure plots don't permanently capture keyboard focus which prevents
  // page scrolling. Remove tabindex attributes and blur non-input focus targets.
  // Also add a document-level keydown handler (added once) to restore
  // page navigation when navigation keys are pressed while a plot has focus.
  let _plotly_scroll_key_handler_added = false;

  function makePlotSafe(plot) {
    try {
      try { console.debug && console.debug('[plotly_scroll] makePlotSafe', plot); } catch (e) {}
      // Remove tabindex so clicks don't trap focus
      if (plot.hasAttribute && plot.hasAttribute('tabindex')) plot.removeAttribute('tabindex');
      // Remove tabindex from descendants that might capture focus
      plot.querySelectorAll && plot.querySelectorAll('[tabindex]').forEach(el => el.removeAttribute('tabindex'));

      // Blur any non-input focus that lands inside the plot
      plot.addEventListener('focusin', (ev) => {
        try {
          try { console.debug && console.debug('[plotly_scroll] focusin inside plot', ev.target); } catch (e) {}
          const t = ev.target;
          if (!t) return;
          if (t.matches && t.matches('input, textarea, [contenteditable]')) return;
          // Defer blur slightly to let the event finish
          setTimeout(() => { try { document.activeElement && document.activeElement.blur && document.activeElement.blur(); } catch (e) {} }, 0);
        } catch (e) {}
      }, true);

      // Add a one-time document keydown handler to restore page navigation
      if (!_plotly_scroll_key_handler_added) {
        _plotly_scroll_key_handler_added = true;
        document.addEventListener('keydown', (e) => {
          try {
            const navKeys = ['PageUp','PageDown','ArrowUp','ArrowDown','Home','End',' '];
            if (navKeys.indexOf(e.key) === -1) return;
            const ae = document.activeElement;
            if (!ae) return;
            if (ae.closest && ae.closest('.plotly-graph-div')) {
              try { ae.blur(); } catch (err) {}
            }
          } catch (err) {}
        }, true);
      }
    } catch (err) {}
  }

  // Attempt to release any pointer captures and blur plot-related focus
  function safeUntrap(plot) {
    try {
      try { console.debug && console.debug('[plotly_scroll] safeUntrap', plot); } catch (e) {}
      // blur focused element if inside plot
      const ae = document.activeElement;
      if (ae && ae.closest && ae.closest('.plotly-graph-div')) {
        try { ae.blur(); } catch (e) {}
      }

      // release any pointer capture on plot or its descendants
      try {
        if (plot && plot.releasePointerCapture && typeof plot.releasePointerCapture === 'function') {
          try { plot.releasePointerCapture && plot.releasePointerCapture(); } catch (e) {}
        }
        plot && plot.querySelectorAll && plot.querySelectorAll('*').forEach(el => {
          try { if (el.releasePointerCapture) el.releasePointerCapture(); } catch (e) {}
        });
      } catch (e) {}

      // dispatch a synthetic mouseup on document to end any drag states
      try {
        const ev = new MouseEvent('mouseup', { bubbles: true, cancelable: true });
        document.dispatchEvent(ev);
        try { console.debug && console.debug('[plotly_scroll] dispatched synthetic mouseup'); } catch (e) {}
      } catch (e) {}
    } catch (err) {}
  }

  // Attach to existing plots
  document.querySelectorAll('.plotly-graph-div').forEach(attachToPlot);
  // Also make them safe for focus/navigation
  document.querySelectorAll('.plotly-graph-div').forEach(makePlotSafe);
  // Ensure we untrap plot interactions on mouseup/mouseleave/modebar clicks
  document.querySelectorAll('.plotly-graph-div').forEach((p) => {
    try {
      p.addEventListener('mouseup', () => safeUntrap(p), { passive: true });
      p.addEventListener('mouseleave', () => safeUntrap(p), { passive: true });
      // modebar buttons have class 'modebar-btn' in many Plotly builds
      p.querySelectorAll && p.querySelectorAll('.modebar, .modebar-btn, .modebar-btn-group').forEach(mb => {
        mb.addEventListener('click', () => safeUntrap(p), { passive: true });
      });
    } catch (e) {}
  });

  // Observe for future plots added dynamically (nbsphinx may inject them)
  const mo = new MutationObserver((records) => {
    for (const r of records) {
      for (const n of r.addedNodes) {
        if (!(n instanceof HTMLElement)) continue;
        if (n.matches && n.matches('.plotly-graph-div')) attachToPlot(n);
        // also check descendants
        n.querySelectorAll && n.querySelectorAll('.plotly-graph-div').forEach((p) => {
          attachToPlot(p); makePlotSafe(p);
          try {
            p.addEventListener('mouseup', () => safeUntrap(p), { passive: true });
            p.addEventListener('mouseleave', () => safeUntrap(p), { passive: true });
            p.querySelectorAll && p.querySelectorAll('.modebar, .modebar-btn, .modebar-btn-group').forEach(mb => {
              mb.addEventListener('click', () => safeUntrap(p), { passive: true });
            });
          } catch (e) {}
        });
      }
    }
  });

  mo.observe(document.body, { childList: true, subtree: true });

  // --- Watchdog: detect stalled page scrolling after plot interaction ---
  // If the user interacts with a plot and a subsequent wheel does not
  // actually scroll the document, call safeUntrap on the last interacted
  // plot to release captures. This is a best-effort mitigation for cases
  // where Plotly prevents default wheel behavior and leaves focus/pointer
  // trapped.
  let lastInteractedPlot = null;

  function markPlotInteraction(plot) {
    lastInteractedPlot = plot;
    try { console.debug && console.debug('[plotly_scroll] markPlotInteraction', plot); } catch (e) {}
    // clear after 2s
    setTimeout(() => { if (lastInteractedPlot === plot) lastInteractedPlot = null; }, 2000);
  }

  // update lastInteractedPlot from various interaction sources
  document.querySelectorAll('.plotly-graph-div').forEach(p => {
    try {
      p.addEventListener('mousedown', () => markPlotInteraction(p), { passive: true });
      p.addEventListener('touchstart', () => markPlotInteraction(p), { passive: true });
      p.addEventListener('focusin', () => markPlotInteraction(p), { passive: true });
    } catch (e) {}
  });

  // Also observe and bind for future plots
  mo.disconnect();
  mo.observe(document.body, { childList: true, subtree: true });
  mo.takeRecords && mo.takeRecords();
  // Keep original mutation logic but ensure new plots get mark handlers
  const mo2 = new MutationObserver((records) => {
    for (const r of records) {
      for (const n of r.addedNodes) {
        if (!(n instanceof HTMLElement)) continue;
        if (n.matches && n.matches('.plotly-graph-div')) {
          attachToPlot(n); makePlotSafe(n);
          try { n.addEventListener('mousedown', () => markPlotInteraction(n), { passive: true }); } catch (e) {}
        }
        n.querySelectorAll && n.querySelectorAll('.plotly-graph-div').forEach((p) => {
          attachToPlot(p); makePlotSafe(p);
          try { p.addEventListener('mousedown', () => markPlotInteraction(p), { passive: true }); } catch (e) {}
        });
      }
    }
  });
  mo2.observe(document.body, { childList: true, subtree: true });

  // Wheel capture-check: compare scrollTop before/after wheel to detect
  // whether the document scrolled. If not, untrap the last interacted plot.
  try {
    document.addEventListener('wheel', (e) => {
      try {
        const doc = document.scrollingElement || document.documentElement || document.body;
        const before = doc.scrollTop;
        setTimeout(() => {
          try {
            const after = doc.scrollTop;
            if (after === before && lastInteractedPlot) {
              try { console.debug && console.debug('[plotly_scroll] wheel did not scroll; calling safeUntrap', lastInteractedPlot); } catch (e) {}
              safeUntrap(lastInteractedPlot);
            }
          } catch (e) {}
        }, 60);
      } catch (e) {}
    }, { passive: true, capture: true });
  } catch (e) {}

})();
