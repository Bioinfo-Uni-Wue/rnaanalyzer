//******************************************************************************
// BED format Checker
//******************************************************************************
var BedChecker = function (handler) {
  "use strict";
  // store a reference to the handler
  this.handler = handler;
  // current parsing function
  this.process = this._process_start;
  // abort flag
  this.give_up = false;
};

BedChecker.prototype._process_start = function (code, type) {
  "use strict";
  return true;
};

// When we're done, call the approprate functions on the handler
BedChecker.prototype._signal_stop = function() {
  if (typeof this.handler.progress == "function") this.handler.progress(1.0);
  if (typeof this.handler.end == "function") this.handler.end();
};

//******************************************************************************
// Public functions
//******************************************************************************

BedChecker.prototype.process_file = function (file_input, sequence_index, prefix) {
  "use strict";
  var me, message, comment_count, line_count, expanded_size, sequence_count, reader;
  me = this;
  comment_count = 0;
  line_count = 0;
  expanded_size = 0;
  sequence_count = 0;
  if (this.give_up) return;

  // file size error
  if (this.handler.max_file_size != null && file_input.size > this.handler.max_file_size) {
    var error = {
      name: "BED File is too large",
      type: FileType.INVALID_BED,
      message: "File size must be no larger than " 
       + (this.handler.max_file_size/1000000) + "MB but is over "
       + (Math.trunc(file_input.size/1000000)) + "MB" + ". "
    }
    this.handler.error_format(error);
    this._signal_stop();
    return;
  }

  if (typeof this.handler.begin == "function") {
    line_count = 0;
    this.handler.begin();
  }

  reader = new FileReader();
  reader.onload = function(evt) {
    "use strict";
    var error = null;
    var enc = new TextEncoder(); 
    if ((error = unusable_format(enc.encode(this.result), 400, file_input.name)) != null) {
      // report error and stop scan as we don't have a chance of understanding this file
      if (typeof me.handler.error_format == "function") 
        me.handler.error_format(error);
      me._signal_stop();
      return;
    }
    // Check for file name too long.
    if (file_input.name.length > me.handler.max_file_name_length) {
      error = {
	name: "Invalid BED file name",
	type: FileType.INVALID_BED,
	message: "The BED file name is too long. "
        + "It should be less than " + me.handler.max_file_name_length
        + " characters but is " + file_input.name.length + "."
      };
    }
    var lines = (this.result).trim().split('\n');
    line_count = lines.length;
    if (typeof me.handler.begin == "function") {
      me.handler.begin(this.size, lines.length);
    }
    for (var line = 0; error == null && line < lines.length; line++) {
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
	    name: "Bad BED file format", 
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
        // check for non-existent sequence name
	if (!sequence_index.has(chrom)) {
	  error = {
	    name: "Bad BED file", 
	    type: FileType.INVALID_BED,
	    message: "Line "  + (line+1) + ": sequence \'" + chrom + "\' not found in genome database."
	  };
	  break;
	}
	var length = Number(sequence_index.get(chrom));
	// The 2nd field should be a non-negative integer
	var start = Number(tabs[1].trim());
        if (isNaN(start)) {
          error = {
            name: "Bad BED file format",
            type: FileType.INVALID_BED,
            message: "Line " + (line+1) + ": field 2 (start) is not a valid number: '" + tabs[1].trim() + "'."
          }
          break;
        }
	if ((Math.floor(start) != start) || start < 0)  {
	  error = {
	    name: "Bad BED file format", 
	    type: FileType.INVALID_BED,
	    message: "Line " + (line+1) + ": field 2 (start) is " + start + ". It must be a non-negative integer."
	  };
	  break;
	}
	// The 3rd field thould be a non-negative integer
	var stop = Number(tabs[2].trim());
        if (isNaN(stop)) {
          error = {
            name: "Bad BED file format",
            type: FileType.INVALID_BED,
            message: "Line " + (line+1) + ": field 3 (stop) is not a valid number: '" + tabs[2].trim() + "'."
          }
          break;
        }
	if ((Math.floor(stop)  != stop) || stop <= 0) {
	  error = {
	    name: "Bad BED file format", 
	    type: FileType.INVALID_BED,
	    message: "Line "  + (line+1) + ": field 3 (stop) is " + stop + ". It must be a non-negative integer."
	  };
          break;
	}
	if (start > stop) {
	  error = {
	    name: "Bad BED file format",
	    type: FileType.INVALID_BED,
	    message: "Line " + (line+1) + ": field 2 (start) is " + start + ". It must be less than field 3 (stop), " + stop + "."
	  };
	  break;
	}
	if (start > length) {
	  error = {
	    name: "Bad BED file", 
	    type: FileType.INVALID_BED,
	    message: "Line " + (line+1) + ": field 2 (start) is " + start + ". It must be less than the chromosome length, " + length + "."
	  };
          break;
	}
	if (stop > length+1) {
	  error = {
	    name: "Bad BED file", 
	    type: FileType.INVALID_BED,
	    message: "Line " + (line+1) + ": field 3 (end) is " + stop + ". It must be no larger than the chromosome length, " + length + "."
	  };
	  break;
	}
        expanded_size += stop - start;
        // Check that expanded BED file has less that max_file_size of sequence.
        if (expanded_size > me.handler.max_file_size) {
	  error = {
	    name: "Bad BED file", 
	    type: FileType.INVALID_BED,
	    message: "The FASTA file created from the BED file will have over " +
             (Math.trunc(expanded_size/1000000)) + "MB but should be no larger than " +
             Math.round(me.handler.max_file_size/1000000) + "MB."
          };
	  break;
        }
	sequence_count++;
      }
    }

    if (lines.length === 0 || comment_count === lines.length) {
      error = {
        name: "Empty BED file",
        type: FileType.INVALID_BED,
        message: "The file doesn't contain any BED data."
      }
    }
    if (error != null && typeof me.handler.error_format == "function") {
      me.handler.error_format(error);
      me._signal_stop();
      return;
    }

    me.handler.line_count = line_count;
    me.handler.comment_count = comment_count;
    me.handler.expanded_size = expanded_size;
    me.handler.sequence_count = sequence_count;
    me._signal_stop();

    // Save expanded size and sequence count in hidden fields for BACK and REFRESH.
    var size_id = prefix + "_expanded_size";
    var count_id = prefix + "_sequence_count";
    document.getElementById(size_id).value = expanded_size;
    document.getElementById(count_id).value = sequence_count;
    me.handler.prefix = prefix;
  };

  if (file_input) {
    reader.readAsText(file_input);
  }
};

BedChecker.prototype.cancel = function () {
  "use strict";
  this.give_up = true;
  this.handler = {};
};

//******************************************************************************
// Bed Handler
//******************************************************************************
var BedHandler = function (options) {
  this.configure(options);
  this.reset();
};

BedHandler.prototype.configure = function (options) {
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

BedHandler.prototype.guess_alphabets = function() {
  "use strict";
  // Only DNA is currently supported for BED files.
  return [AlphStd.DNA];
};

BedHandler.prototype.reset = function () {
  // have the file details changed?
  this.updated = false;
  // the part of the file processed
  this.fraction = 0;
  // BED details
  this.file_size = 0;
  this.line_count = 0;
  this.comment_count = 0;
  this.expanded_size = 0;
  this.sequence_count = 0;
  // keep track of problems found
  this.error_type = 0;
  this.error_name = null;
  this.error_message = null;
  this.encoding_error = null;
  //this.missing_name = new FileFaults();
};

BedHandler.prototype.summary = function () {
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
  // file size warning
  if (this.max_file_size != null && this.file_size > this.max_file_size) {
    add(false, "FASTA File is too large. ", ["File is over " + (Math.trunc(this.file_size/1000000)) + "MB but should be no larger than " + 
      (this.max_file_size/1000000) + "MB."]);
  }
  help = " - re-save as plain text; either Unicode UTF-8 (no Byte Order Mark) or ASCII";
  if (this.error_name != null) {
    switch (this.error_type) {
      case FileType.ENCODING:
        add(true, "Bad encoding \"" + this.error_name + "\"" + help);
        break;
      case FileType.BINARY:
        add(true, "Bad format \"" + this.error_name + "\"" + help);
        break;
      case FileType.COMPRESSED:
        add(true, "Bad format \"" + this.error_name + "\" - must be decompressed first");
        break;
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
BedHandler.prototype.progress = function (fraction) {
  "use strict";
  this.fraction = fraction;
};

// Reading of the file has begun
BedHandler.prototype.begin = function (num_lines, file_size) {
  "use strict";
  this.reset();
  this.line_count = num_lines;
  this.file_size = file_size;
  this.updated = true;
};

// Reading of the file has finished (perhaps early due to an error)
BedHandler.prototype.end = function () {
  "use strict";
  this.updated = true;
};

// Parsing has stopped due to an unreadable file format
BedHandler.prototype.error_format = function (error) {
  "use strict";
  this.error_type = error.type;
  this.error_name = error.name;
  this.error_message = error.message
  this.updated = true;
};

