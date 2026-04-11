/* ==========================================================================
   adataviz – custom enhancements for the Furo-based documentation
   ========================================================================== */

(function () {
  "use strict";

  /* --- Back-to-top button ------------------------------------------------ */
  var btn = document.createElement("button");
  btn.id = "back-to-top";
  btn.setAttribute("aria-label", "Back to top");
  btn.innerHTML = "&#8679;"; // ⇧
  document.body.appendChild(btn);

  var scrollTarget =
    document.querySelector(".main") || document.documentElement;

  function toggleBtn() {
    var scrollY =
      scrollTarget.scrollTop || document.documentElement.scrollTop || 0;
    if (scrollY > 400) {
      btn.classList.add("visible");
    } else {
      btn.classList.remove("visible");
    }
  }

  // Furo scrolls .main rather than the window on wide viewports
  scrollTarget.addEventListener("scroll", toggleBtn, { passive: true });
  window.addEventListener("scroll", toggleBtn, { passive: true });

  btn.addEventListener("click", function () {
    scrollTarget.scrollTo({ top: 0, behavior: "smooth" });
    window.scrollTo({ top: 0, behavior: "smooth" });
  });

  /* --- External links open in new tab ------------------------------------ */
  document.querySelectorAll('a[href^="http"]').forEach(function (a) {
    if (!a.hostname || a.hostname === window.location.hostname) return;
    a.setAttribute("target", "_blank");
    a.setAttribute("rel", "noopener noreferrer");
  });

  /* --- Wrap tables in a scrollable container ----------------------------- */
  document.querySelectorAll(".content table").forEach(function (table) {
    if (table.parentElement.classList.contains("table-wrapper")) return;
    var wrapper = document.createElement("div");
    wrapper.className = "table-wrapper";
    wrapper.style.overflowX = "auto";
    table.parentNode.insertBefore(wrapper, table);
    wrapper.appendChild(table);
  });
})();
