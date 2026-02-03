// Prevent page scroll when mouse is over a Plotly graph and forward wheel events
// to the inner SVG so Plotly can perform scroll-to-zoom.

(function () {
  function onWheel(e) {
    // Find a closest Plotly graph container
    const plot = e.target.closest('.plotly-graph-div');
    if (!plot) return;

    // Prevent page from scrolling
    e.preventDefault();

    // Try to forward the wheel event to the inner SVG (where Plotly listens)
    const svg = plot.querySelector('svg');
    if (svg) {
      const evt = new WheelEvent('wheel', {
        deltaX: e.deltaX,
        deltaY: e.deltaY,
        clientX: e.clientX,
        clientY: e.clientY,
        bubbles: true,
        cancelable: true,
        composed: true
      });
      svg.dispatchEvent(evt);
    }
  }

  // Use capture and passive:false so we can call preventDefault()
  document.addEventListener('wheel', onWheel, { passive: false, capture: true });
})();
