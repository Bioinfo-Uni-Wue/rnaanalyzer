/******************************************************************************
 * Create a LOCUS file upload indicator handler.
 ******************************************************************************/
var LOCUSIH = function (container, indicator, popup, options, change_handler) {
  "use strict";
  this.container = container;
  this.indicator = indicator;
  this.popup = popup;
  this.change_handler = change_handler;
  this.update_timer = null;
  this.display_type = 0;
  this.display_error = true;
  this.display_warn = false;
  LociHandler.call(this, options);
};

LOCUSIH.prototype = Object.create(LociHandler.prototype);
LOCUSIH.prototype.constructor = LOCUSIH;

/******************************************************************************
 * Clears all the indicators.
 ******************************************************************************/
LOCUSIH.prototype.clear_indicators = function() {
  "use strict";
  if (this.update_timer) {
    clearInterval(this.update_timer);
    this.update_timer = null;
  }
  this.popup.innerHTML = "";
  this.indicator.style.width = "0%";
  substitute_classes(this.container, ["good", "error", "warning"], []);
  this.display_type = 0;
};

LOCUSIH.prototype._update_display = function () {
  "use strict";
  var table, summary, i, j, msg;
  var row, cell, i, ul, li;

  this.indicator.style.width = (this.fraction * 100) + "%";
  if (!this.updated) return;

  summary = this.summary();
  if (summary.error) {
    if (this.display_type != 3) {
      substitute_classes(this.container, ["good", "warning"], ["error"]);
      this.display_type = 3;
    }
  } else if (summary.warning) {
    if (this.display_type != 2) {
      substitute_classes(this.container, ["good", "error"], ["warning"]);
      this.display_type = 2;
    }
  } else {
    if (this.display_type != 1) {
      substitute_classes(this.container, ["warning", "error"], ["good"]);
      this.display_type = 1;
    }
  }
  table = document.createElement("table");
  for (i = 0; i < summary.messages.length; i++) {
    msg = summary.messages[i];
    row = table.insertRow(table.rows.length);
    row.className = (msg.is_error ? "error" : "warning");
    cell = row.insertCell(row.cells.length);
    cell.appendChild(document.createTextNode(msg.is_error ? "\u2718" : "\u26A0"));
    cell = row.insertCell(row.cells.length);
    cell.appendChild(document.createTextNode(msg.message));
    if (typeof msg.reasons !== "undefined") {
      ul = document.createElement("ul");
      for (j = 0; j < msg.reasons.length; j++) {
        li = document.createElement("li");
        li.appendChild(document.createTextNode(msg.reasons[j]));
        ul.appendChild(li);
      }
      cell.appendChild(ul);
    }
  }

  this.popup.innerHTML = "";
  this.popup.appendChild(table);

}

LOCUSIH.prototype.reset = function () {
  "use strict";
  LociHandler.prototype.reset.apply(this, arguments);
  this.clear_indicators();
};

LOCUSIH.prototype.begin = function () {
  "use strict";
  var me;
  LociHandler.prototype.begin.apply(this, arguments);
  if (this.update_timer) clearInterval(this.update_timer);
  me = this; // reference 'this' inside closure
  this.update_timer = setInterval(function () { me._update_display();}, 100);
  this.change_handler();
};

/******************************************************************************
 *
 ******************************************************************************/
LOCUSIH.prototype.end = function () {
  "use strict";
  LociHandler.prototype.end.apply(this, arguments);
  if (this.update_timer) {
    clearInterval(this.update_timer);
    this.update_timer = null;
  }
  this._update_display();
  if (this.display_type == 1) this.clear_indicators();
  this.change_handler();
};

/******************************************************************************
 * Construct the manager of the loci input.
 ******************************************************************************/
var LociInput = function(container, options) {
  "use strict";
  var me, box;
  // make 'this' accessable in inner scopes
  me = this;
  // store the parameters
  this.container = container;
  this.options = options;
  // lookup relevent components
  // get the locus file related parts
  this.file_surround = this.container.querySelector("span.loci_file");
  this.file_indicator = this.file_surround.querySelector("span.indicator");
  this.file_input = this.file_surround.querySelector("input");
  this.file_popup = this.file_surround.querySelector("div.popup");
  this.file_dbh = new LOCUSIH(this.file_surround, this.file_indicator, this.file_popup, options, function() {
      me._fire_loci_checked_event();
    });
  // create a list of submittable fields so we can disable the ones we are not using
  this.submittables = [this.file_input];
  // other things
  this.parser = null;
  this.timer = null;
  // initialise
  this._file_update();
  // add listeners
  // detect file changes
  this.file_input.addEventListener('change', function() {
    me._file_update();
  }, false);
};

/******************************************************************************
 *
 ******************************************************************************/
LociInput.prototype.check = function() {
  "use strict";
  var exists, summary;
  exists = this.file_input.value.length > 0;
  if (!exists) {
    alert("Please input the " + this.options.field + ".");
    return false;
  }
  summary = this.file_dbh.summary();
  if (summary.error) {
    alert("Please correct errors in the " + this.options.field + ".");
    return false;
  }
  return true;
};

/******************************************************************************
 *
 ******************************************************************************/
// reset the fields to a consistent state
LociInput.prototype.reset = function() {
  "use strict";
  if (this.parser) this.parser.cancel();
  this.parser = null;
  if (this.timer) window.clearTimeout(this.timer);
  this.timer = null;
  this.file_input.value = "";
  this.file_dbh.reset();
};

/******************************************************************************
 *
 ******************************************************************************/
LociInput.prototype._fire_loci_checked_event = function() {
  "use strict";
  var me;
  me = this;
  try {
    // IE sometimes has problems with this line.
    // I think they are related to the page not being fully loaded.
    this.container.dispatchEvent(new CustomEvent("loci_checked", {detail: {controler: me}}));
  } catch (e) {
    if (e.message && e.name && window.console) {
      console.log("Suppressed exception " + e.name + ": " + e.message);
    }
  }
};

/******************************************************************************
 *
 ******************************************************************************/
LociInput.prototype._file_update = function() {
  "use strict";
   this._file_validate();
};

/******************************************************************************
 *
 ******************************************************************************/
LociInput.prototype._file_validate = function() {
  var file;
  if (this.parser) this.parser.cancel();
  //file = this.file_input.files[0];
  if (!(file = this.file_input.files[0])) {
    this.file_dbh.clear_indicators();
    return;
  }
  this.file_dbh.configure(this.options);
  this.parser = new LociChecker(this.file_dbh);
  this.parser.process_file(file);
};
