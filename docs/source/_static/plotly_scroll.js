// Prevent page scroll when mouse is over a Plotly graph and forward wheel events
// to the inner SVG so Plotly can perform scroll-to-zoom.

(function () {
  function onWheel(e) {
    // Find a closest Plotly graph container
    const plot = e.target.closest('.plotly-graph-div');
    if (!plot) return;
    // Ignore events that were forwarded by this script to avoid recursion
    if (e.__plotly_forwarded) return;

    // Prevent page from scrolling
    e.preventDefault();

    // Try to forward the wheel event to the inner SVG (where Plotly listens)
    const svg = plot.querySelector('svg');
    if (svg) {
      // Create a non-bubbling event so it won't be captured again by our
      // document-level listener (which would cause infinite recursion).
      const evt = new WheelEvent('wheel', {
        deltaX: e.deltaX,
        deltaY: e.deltaY,
        deltaMode: e.deltaMode || 0,
        clientX: e.clientX,
        clientY: e.clientY,
        bubbles: false,
        cancelable: true,
        composed: false
      });

      // Mark the synthetic event so our handler ignores it if seen.
      try {
        Object.defineProperty(evt, '__plotly_forwarded', { value: true });
      } catch (err) {
        // Some browsers may not allow defining custom props — as a fallback
        evt.__plotly_forwarded = true;
      }

      // Give focus to the graph so Plotly receives scroll-zoom input
      try { svg.focus && svg.focus(); } catch (err) {}

      svg.dispatchEvent(evt);
    }
  }

  // Use capture and passive:false so we can call preventDefault()
  document.addEventListener('wheel', onWheel, { passive: false, capture: true });
})();
