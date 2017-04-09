# prcnToBed v1.7, 3 Feb 2015
# Convert PERCON DISTANCE output file to BED (Browser Extensible Display) format.
# Author: Lev I. Uralsky (Institute of Molecular Genetics, Moscow, Russia)
# 
# v1.0, 18 Nov 2012 - initial release
# v1.1, 23 Nov 2012
# v1.2, 30 Nov 2012
# v1.3,  1 Mar 2013 - adding cromosome ids conversion
# v1.4, 15 Mar 2014 - added new colors
# v1.5, 21 Oct 2014 - changes to get chromosome name
# v1.6, 12 Jan 2015 - corrected syntax; added commentaries;
#                     removing unused elements
# v1.7,  3 Feb 2015 - improvements for better handling different
#                     versions of the outputs
#
# Usage: gawk -f prcnToBed.awk malux3 > out.bed

BEGIN {
  # ignoring the letter case in monomers types
  IGNORECASE = 1;
  # setting output field separator to tab
  OFS = "\t";
  # constant with file name, for converting ids into the chromosome names
  idsFile = "id.txt";
  # initializing variable used for flagging in converting ids
  foundCM = 0;
}

# this is the first line in the report output
/            PERCON/ {
  # sets the variable identifier
  seqID = "HEADER";
}

# zero out a 'seqID' variable in these lines
/^seq ID/ ||
/subtypes distribution/ ||
/SF unidentified/ ||
/^         type/ {
  # saves a previous value
  if (seqID) {
    prevSeqID = seqID;
  }
  seqID = "";
  # moves on to the next line
  getline;
}

# this is the line that begins footer in the report output
/invalid format doks/ ||
/total      calculation time/ {
  # sets the variable identifier
  seqID = "FOOTER";
}

# at this point line must begin with sequence accession or chromosome ids
!seqID && /^[^[:space:]]/ && !/^>/ && /[^ACGTN]/ {
  # skipping if the identifier is the same as before
  if ($1 != prevSeqID) {
    # if not save it
    seqID = $1;

    # deleting character '>' if it present
    sub(/\76/, "", seqID);

    # if field begin with chromosome id
    if (seqID ~ /^[1-9XY]/) {
      # adding "chr" string to the beginning
      chrid = "chr" seqID;
    }

    # if field beginning with special id
    if (seqID ~ /^CM|^FR|^NW_/) {
      # if we are sure what this operation performed only once
      if (!(foundCM)) {
        # and file with ids exists
        if (exists(idsFile)) {
          # resetting counter
          rec = 0;
          # set __input variable from next record of idsFile
          while ((getline __input < idsFile) > 0) {
            # adding string from idsFile to array indexed by line number
            aCM[++rec] = __input;
          }
          # closing the file after adding the last line
          close(idsFile);
          # continuing with results as array
          for (i in aCM) {
            # splitting each item by tab into two part
            split(aCM[i], aCMLine, "\t");
            # first one contain id into which you want to convert
            chr = aCMLine[1];
            # second one contain id that you want to convert
            cm  = aCMLine[2];
            # create new array indexed by second id
            aChrCM[cm] = chr;
          }
        } else {
          # if file for converting does not exist, prints error message
          # but continue executing
          print "Error! File " idsFile " for conversion CM id's does not exist...";
        }
        # set flag to shows that the operation has already been accomplished once
        foundCM = 1;
      } 
    }

    # sets id from current line
    currCM = seqID;
    # and search in array for converting
    if (currCM in aChrCM) {
      # if it found sets the variable which beginning with string "chr"
      # and chromosome id
      chrid = "chr" aChrCM[currCM];
    } else {
      # otherwise sets the variable to current id, as is
      chrid = seqID;
    }
  }
}

# now we processing line, which contain information about the monomer
/^         \w/ &&
!/type/ &&
seqID != "FOOTER" {
  # default setting for strand variable is plus
  strand = "+";
  # sets variable for monomer type
  type1 = $1;
  # sets variable for box type in upper case
  box = toupper($2);

  # Here we deal with the complex case, which occurs when the coordinates are
  # represented by eight-digit value, and two fields stick together in one.
  # in the usual case field length are less than 8 symbols
  if (length($3) < 10) {
    # and we set a variable with the initial and final coordinates of the monomer
    begin = $3;
    end = $4;
  } else {
    # but in the complex case, first we need to find the middle of string
    mid3 = length($3)/2;
    # then initial coordinate will be in the first part of the string
    begin = substr($3, 1, mid3);
    # and final coordinate at the last part
    end = substr($3, mid3 + 1);
  }

  # in case, the line contains the letter C, which shows, what strand is reverse
  if ($0 ~ /C[^a-z]/) {
    # sets the strand variable to minus
    strand = "-";
  }

  # handles special cases with M1 subfamilies
  if (length($NF) == 2) {
    # sets the variable, if the name of the subfamily M1 is
    # present in the last column
    mtype = $NF;
  } else if ((length($(NF-1)) == 2) && ($(NF-1) != "..")) {
    # or if the name of the subfamily M1 is present in the last
    # column before and is not an '..'
    mtype = $(NF-1);
  } else {
    # if the name of the subfamily M1 is not found
    mtype = "";
  }

  # handles special cases with SF classes (letter case ignored)
  if ((type1 ~ /M1/) && (mtype ~ /Ua/)) {
    # does not assign variable in M1 monomers with subtype Ua
    type1 = "";
  } else if ((type1 ~ /M1/) && (mtype != "")) {
    # does not assign variable in M1 monomers
    type1 = "";
  } else if (mtype ~ /Ua/) {
    # does not assign variable in Ua monomers
    mtype = "";
  }

  # here is assigned a color for each type
  switch (type1 "" mtype) {
    case /M1/:
      color = "255,255,0";
      break;
    case /Aa/:
      color = "172,172,172";
      break;
    case /Ja/:
      color = "225,126,231";
      break;
    case /Ba/:
      color = "255,153,0";
      break;
    case /Ca/:
      color = "224,0,64";
      break;
    case /Ia/:
      color = "153,51,102";
      break;
    case /Oa/:
      color = "255,229,153";
      break;
    case /Na/:
      color = "32,160,64";
      break;
    case /Ka/:
      color = "0,255,0";
      break;
    case /Ea/:
      color = "153,153,255";
      break;
    case /Ha/:
      color = "127,96,0";
      break;
    case /Fa/:
      color = "0,128,128";
      break;
    case /Ga/:
      color = "255,255,0";
      break;
    case /R1/:
      color = "0,96,192";
      break;
    case /R2/:
      color = "102,153,255";
      break;
    case /D1/:
      color = "153,0,255";
      break;
    case /D2/:
      color = "210,110,250";
      break;
    case /W/:
      color = "0,255,255";
      break;
    case /J1/:
      color = "234,153,153";
      break;
    case /J2/:
      color = "255,204,204";
      break;
    default:
      color = "85,85,85";
  }

  # printing the resulting line for each monomer
  print chrid, begin, end, type1 mtype "," box, "0", strand, begin, end, color;

}

# function which checks existence of file
function exists(file,  __tmp) {
  if ((getline __tmp < file) > 0) {
    close(file);
    return(1);
  }
  return(0);
}
