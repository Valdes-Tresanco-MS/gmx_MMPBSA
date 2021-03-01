/*
 *                          GPLv3 LICENSE INFO
 *
 * Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco
 * Copyright (C) 2014  Jason Swails, Bill Miller III, and Dwight McGee
 *
 *  Project: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA
 *
 *  This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 */

var topbtn = document.getElementById("topbtn");
var footer = document.getElementsByTagName('footer')

function isScrolledIntoView(element) {
    var elemRect = element.getBoundingClientRect();
    return (elemRect.top >= 0) && (elemRect.bottom <= window.innerHeight);
}

window.onscroll = function() {scrollFunction()};

function scrollFunction() {
  if ((document.body.scrollTop > 30 || document.documentElement.scrollTop > 30) && !(isScrolledIntoView(footer[0])) ) {
    topbtn.style.display = "block";
  } else {
    topbtn.style.display = "none";
  }
}

// Scroll to top
function topFunction() {
  document.body.scrollTop = 0;
  document.documentElement.scrollTop = 0;
}