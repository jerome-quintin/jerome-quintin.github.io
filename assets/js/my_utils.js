"use strict";

// Based on https://stackoverflow.com/a/22480938/1695428
function isCompletelyScrolledIntoView(el) {
  var rect = el.getBoundingClientRect();
  var elemTop = rect.top;
  var elemBottom = rect.bottom;

  // Only completely visible elements return true:
  var isVisible = (elemTop >= 0) && (elemBottom <= window.innerHeight);

  return isVisible;
}

function isPartiallyScrolledIntoView(el) {
  var rect = el.getBoundingClientRect();
  var elemTop = rect.top;
  var elemBottom = rect.bottom;

  // Partially visible elements return true:
  var isVisible = (elemTop < window.innerHeight) && elemBottom >= 0;

  return isVisible;
}
