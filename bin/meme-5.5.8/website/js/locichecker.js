var LociCheckerUtil = {};

//******************************************************************************
// Loci (BED) format Checker 
//******************************************************************************
var LociChecker = function (handler) {
  "use strict";
  // store a reference to the handler
  this.handler = handler;
  // current parsing function
  this.process = this._process_start;
  // abort flag
  this.give_up = false;
};

LociChecker.prototype._process_start = function (code, type) {
  "use strict";
  return true;
};

// When we're done, call the approprate functions on the handler
LociChecker.prototype._signal_stop = function() {
  if (typeof this.handler.progress == "function") this.handler.progress(1.0);
  if (typeof this.handler.end == "function") this.handler.end();
};

//******************************************************************************
// Public functions
//******************************************************************************

LociChecker.prototype.process_file = function(file_input) {
  "use strict";
  var me, message, comment_count, reader;
  me = this;
  var error = {};
  comment_count = 0;
  if (this.give_up) return;

  // file size error 
  if (this.handler.max_file_size != null && file_input.size > this.handler.max_file_size) {
    error = {
      name: "Loci BED File is too large",
      type: FileType.INVALID_BED,
      message: "File size must be no larger than " + (this.handler.max_file_size/1000000) + "MB but is over " 
       + (Math.trunc(file_input.size/1000000)) + "MB" + ". "
    }
    this.handler.error_format(error);
    this._signal_stop();
    return;
  }

  if (typeof this.handler.begin == "function")
    this.handler.begin(file_input.name, file_input.size);

  reader = new FileReader();
  reader.onload = function(evt) {
    "use strict";
    // Check for BED file name too long.
    if (file_input.name.length > me.handler.max_file_name_length) {
      error = {
        name: "Bad Loci BED file name",
        type: FileType.INVALID_BED,
        message: "The Loci BED file name is too long. "
        + "It should be less than " + me.handler.max_file_name_length
	+ " characters but is " + file_input.name.length + "."
      };
    }
    var lines = (this.result).trim().split('\n');
    for (var line = 0; error !== {} && line < lines.length; line++) {
      if (typeof me.handler.progress == "function") {
        me.handler.progress(line / lines.length);
      }
      // Skip blank lines, comment lines, browser lines and track lines.
      if (lines[line].trimStart() == "" ||
	lines[line].charAt(0) == '#' ||
	lines[line].startsWith("browser") ||
	lines[line].startsWith("track")) {
        comment_count++;
      } else {
	var tabs = lines[line].split('\t');
	if (tabs.length < 3 || tabs.length > 12) {
	  error = {
	    name: "Bad Loci BED file format",
	    type: FileType.INVALID_BED,
            message: "Invalid format found on line " + (line + 1) + ". "
	      + "It has " + tabs.length + " tab delimited fields. "
	      + "It should have from 3 to 12 fields."
          };
          break;
	}
        var chrom = tabs[0].trim();
        // check for too long sequence name
        if (me.handler.max_name_length != null && chrom.length > me.handler.max_name_length) {
          error = {
            name: "Bad BED file format",
            type: FileType.INVALID_BED,
            message: "Line " + (line+1) + ": sequence name is too long. It should be less than " +
              me.handler.max_name_length + " characters but is " + chrom.length + "."
          };
          break;
        }
	// The 2nd field must be a non-negative integer.
	var start = Number(tabs[1].trim());
	if (isNaN(start)) {
	  error = {
	    name: "Bad BED file format",
	    type: FileType.INVALID_BED,
	    message: "Line " + (line+1) + ": field 2 (start) is not a valid number: '" + tabs[1].trim() + "'."
          }
	  break;
        }
	if ((Math.floor(start) != start) || start < 0) {
	  error = {
	    name: "Bad BED file format",
	    type: FileType.INVALID_BED,
	    message: "Line " + (line+1) + ": field 2 (start) is " + start + ". It must be a non-negative integer."
          }
	  break;
	}
	// The 3rd field must be a non-negative integer.
	var stop = Number(tabs[2].trim());
	if (isNaN(stop)) {
	  error = {
	    name: "Bad BED file format",
	    type: FileType.INVALID_BED,
	    message: "Line " + (line+1) + ": field 3 (stop) is not a valid number: '" + tabs[2].trim() + "'."
          }
	  break;
        }
	if ((Math.floor(stop) != stop) || stop < 0) {
	  error = {
	    name: "Bad BED file format",
	    type: FileType.INVALID_BED,
	    message: "Line " + (line+1) + ": field 3 (stop) is " + stop + ". It must be a non-negative integer."
          }
	  break;
	}
	if (start > stop) {
	  error = {
	    name: "Bad Loci BED file format",
	    type: FileType.INVALID_BED,
	    message: "Line " + (line+1) + ": field 2 (start) is " + start + ". It must be less than field 3 (stop), " + stop + "."
          }
	  break;
	}
      }
    }

    if (lines.length === 0 || comment_count === lines.length) {
      error = {
	name: "Empty Loci BED file",
	type: FileType.INVALID_BED,
        message: "The file doesn't contain any BED data."
      }
    }
    if (error !== {} && typeof me.handler.error_format == "function") {
      me.handler.error_format(error);
      me._signal_stop();
      return;
    }
  };

  if (file_input) {
    reader.readAsText(file_input);
  }
};

LociChecker.prototype.cancel = function () {
  "use strict";
  this.give_up = true;
  this.handler = {};
};

//******************************************************************************
// Loci Handler
//******************************************************************************
var LociHandler = function (options) {
  this.configure(options);
  this.reset();
};

LociHandler.prototype.configure = function (options) {
  "use strict";
  if (typeof options != "object" || options == null) options = {};
  // specify a maximum file size
  if (typeof options.max_file_size == "number" && options.max_file_size >= 1) {
    this.max_file_size = options.max_file_size;
  } else {
    this.max_file_size = null;
  }
  // specify a maximum sequence name length
  if (typeof options.max_name_len == "number" && options.max_name_len >= 1) {
    this.max_name_length = options.max_name_len;
  } else {
   this.max_name_length = null;
  }
  // specify a maximum file name length
  if (typeof options.max_file_name_len == "number" && options.max_file_name_len >= 1) {
    this.max_file_name_length = options.max_file_name_len;
  } else {
    this.max_file_name_length = null;
  }
};

LociHandler.prototype.reset = function () {
  // have the file details changed?
  this.updated = false;
  // the part of the file processed
  this.fraction = 0;
  // fasta details
  this.file_size = 0;
  this.file_symbols = "";
};

LociHandler.prototype.summary = function () {
  "use strict";
  var error, warning, messages, reason, reasons, letters, add;
  var help;
  // setup
  error = false;
  warning = false;
  messages = [];
  // create closure to add messages
  add = function(is_error, message, reasons) {
    "use strict";
    messages.push({"is_error": is_error, "message": message, "reasons": reasons});
    if (is_error) error = true;
    else warning = true;
  };
  if (this.error_name != null) {
    switch(this.error_type) {
      case FileType.INVALID_BED:
        add(true, this.error_name + ": " + this.error_message);
        break;
    }
  }
  // clear updated state
  this.updated = false;
  // return state
  return {"error": error, "warning": warning, "messages": messages};
};

// tracks the progress of reading the file
LociHandler.prototype.progress = function (fraction) {
  "use strict";
  this.fraction = fraction;
};

// Reading of the file has begun
LociHandler.prototype.begin = function (file_name, file_size) {
  "use strict";
  this.reset();
  this.file_name = file_name;
  this.file_size = file_size;
  this.updated = true;
};

// Reading of the file has finished (perhaps early due to an error)
LociHandler.prototype.end = function () {
  "use strict";
  this.updated = true;
};

// Parsing has stopped due to an unreadable file format
LociHandler.prototype.error_format = function (error) {
  "use strict";
  this.error_type = error.type;
  this.error_name = error.name;
  this.error_message = error.message
  this.updated = true;
};
