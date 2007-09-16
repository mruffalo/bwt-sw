#!/usr/bin/env perl

# Location of binaries
$blastloc = "/data1/blast/bin/";
$bwtswloc = "./";
$compareloc = "./";

# Location of databases
$blastdbloc = "/data1/blast/bin/";
$bwtswdbloc = "./";

# Name of binaries
$blast = "blastall";
$bwtsw = "bwtsw";
$compare = "blastcompare";

# Extensions
$bwtswResult = ".bwtsw";
$blastResult = ".blast";
$blastMissed = ".blast.missed";
$bwtswMissed = ".bwtsw.missed";
$rankByQuery = ".rank.query";
$evalueByQuery = ".evalue.query";
$rankSummary = ".rank.summary";
$evalueSummary = ".evalue.summary";
$align = ".align";

# Temporary files
$tempbwtswResult = ".t.bwtsw";
$tempblastResult = ".t.blast";
$tempbwtswMissed = ".t.bwtsw.m";
$tempblastMissed = ".t.blast.m";

# Comparison parameter
$compareRankParameter = "-r 100";
$compareEvalueParameter = "-e 1e-100";

# check arguments
if ($#ARGV < 2) {
  print $0 . " <database name> <query name> <output name> <arguments>\n";
  exit;
}
for ($i = 3; $i <= $#ARGV; $i++) {
  $arg = $arg . " " . $ARGV[$i];
}

# Incorporate location settings
if (substr($blastloc, -1) eq "/") {
  $blast = $blastloc . $blast;
} else {
  $blast = $blastloc . "/" . $blast;
}
if (substr($bwtswloc, -1) eq "/") {
  $bwtsw = $bwtswloc . $bwtsw;
} else {
  $bwtsw = $bwtswloc . "/" . $bwtsw;
}
if (substr($compareloc, -1) eq "/") {
  $compare = $compareloc . $compare;
} else {
  $compare = $compareloc . "/" . $compare;
}
if (substr($blastdbloc, -1) eq "/") {
  $blastdb = $blastdbloc . $ARGV[0];
} else {
  $blastdb = $blastdbloc . "/" . $ARGV[0];
}
if ($bwtswdbloc !~ /\S/) {
  $bwtswdb = $ARGV[0];
} else {
  if (substr($bwtswdbloc, -1) eq "/") {
    $bwtswdb = $bwtswdbloc . $ARGV[0];
  } else {
    $bwtswdb = $bwtswdbloc . "/" . $ARGV[0];
  }
}

# Execute bwtsw
$command = $bwtsw . " -m 9 " . $bwtswdb . " " . $ARGV[1] . " " . $ARGV[2] . $bwtswResult . " -align " . $ARGV[2] . $align . $arg;
print $command . "\n";
$rc = system($command);
if ($rc != 0) {
  die("bwtsw failed. Script terminated.");
}

# Execute blast
$command = $blast . " -p blastn -m 9 -d " . $blastdb . " -i " . $ARGV[1] . " -o " . $ARGV[2] . $blastResult . $arg;
print $command . "\n";
$rc = system($command);
if ($rc != 0) {
  die("blast failed. Script terminated.");
}

print "Comparing results..";

# Open files
open(QUERY, $ARGV[1]) or die("Cannot open query!\n");
open(BLAST, $ARGV[2] . $blastResult) or die("Cannot open blast result!\n");
open(BWTSW, $ARGV[2] . $bwtswResult) or die("Cannot open bwtsw result!\n");
open(BLAST_MISSED, ">" . $ARGV[2] . $blastMissed) or die("Cannot open blast missed alignment!\n");
open(BWTSW_MISSED, ">" . $ARGV[2] . $bwtswMissed) or die("Cannot open bwtsw missed alignment!\n");

# Delete summary file
unlink($ARGV[2] . $rankByQuery);
unlink($ARGV[2] . $evalueByQuery);

# Loop through query file
$c = 1;
while ($queryName = <QUERY>) {
  if (substr($queryName, 0, 1) eq ">") {
    
    $queryName = substr(substr($queryName, 1), 0, -1);
    
    # open temporary files
    open(TEMP_BLAST_RESULT, ">" . $tempblastResult) or die("Cannot open temp blast result!");
    open(TEMP_BWTSW_RESULT, ">" . $tempbwtswResult) or die("Cannot open temp bwtsw result!");;

    # Retrieve blast result
    $blast_c = $blast_t;
    while ($blast_t = <BLAST> and substr($blast_t, 0, 1) eq "#") {
      $blast_c = $blast_c . $blast_t;
      if (substr($blast_t, 0, 8) eq "# Query:") {
        $blastQueryName = substr(substr($blast_t, 9), 0, -1);
      }
    }
    if ($blastQueryName eq $queryName) {
      print TEMP_BLAST_RESULT $blast_t;
      while ($blast_t = <BLAST> and substr($blast_t, 0, 1) ne "#") {
        print TEMP_BLAST_RESULT $blast_t;
      }
    }
    
    # Retrieve bwtsw result
    $bwtsw_c = $bwtsw_t;
    while ($bwtsw_t = <BWTSW> and substr($bwtsw_t, 0, 1) eq "#") {
      $bwtsw_c = $bwtsw_c . $bwtsw_t;
    }
    print TEMP_BWTSW_RESULT $bwtsw_t;
    while ($bwtsw_t = <BWTSW> and substr($bwtsw_t, 0, 1) ne "#") {
      print TEMP_BWTSW_RESULT $bwtsw_t;
    }

    close(TEMP_BLAST_RESULT);
    close(TEMP_BWTSW_RESULT);
    
    # Compare result by evalue
    $command = $compare . " " . $compareEvalueParameter . " " . $tempblastResult . " " . $tempbwtswResult . " " . "\"" . $queryName . "\"" . " " . $tempbwtswMissed . " " . $tempblastMissed . " >> " . $ARGV[2] . $evalueByQuery;
    print $command . "\n";
    $rc = system($command);
    if ($rc != 0) {
      die("blastcompare failed. Script terminated.");
    }
    # Compare result by rank
    $command = $compare . " " . $compareRankParameter . " " . $tempblastResult . " " . $tempbwtswResult . " " . "\"" . $queryName . "\"" . " " . $tempbwtswMissed . " " . $tempblastMissed . " >> " . $ARGV[2] . $rankByQuery;
    print $command . "\n";
    $rc = system($command);
    if ($rc != 0) {
      die("blastcompare failed. Script terminated.");
    }

    # Check output of comparison
    open(TEMP_BWTSW_MISSED, $tempbwtswMissed) or die("Cannot open temp bwtsw missed!");
    $t = <TEMP_BWTSW_MISSED>;
    if ($t =~ /\S/) {
       print BWTSW_MISSED $blast_c . $t;
       while ($t = <TEMP_BWTSW_MISSED>) {
          print BWTSW_MISSED $t;
       }
    }
    close(TEMP_BWTSW_MISSED);
    
    print $tempblastMissed;

    open(TEMP_BLAST_MISSED, $tempblastMissed) or die("Cannot open temp blast missed!");
    $t = <TEMP_BLAST_MISSED>;
    if ($t =~ /\S/) {
        print BLAST_MISSED $bwtsw_c . $t;
        while ($t = <TEMP_BLAST_MISSED>) {
            print BLAST_MISSED $t; 
        }
    }
    close(TEMP_BLAST_MISSED);

    $c = $c + 1;
    
  }
}

# Close files
close(QUERY);
close(BLAST);
close(BWTSW);
close(BLAST_MISSED);
close(BWTSW_MISSED);

# Delete temp files
unlink($tempbwtswResult);
unlink($tempblastResult);
unlink($tempbwtswMissed);
unlink($tempblastMissed);

# Summarize
$command = $compare . " -s " . $compareEvalueParameter . " < " . $ARGV[2] . $evalueByQuery . " > " . $ARGV[2] . $evalueSummary;
print $command . "\n";
$rc = system($command);
if ($rc != 0) {
    die("blastcompare failed. Script terminated.");
}

$command = $compare . " -s " . $compareRankParameter . " < " . $ARGV[2] . $rankByQuery . " > " . $ARGV[2] . $rankSummary;
print $command . "\n";
$rc = system($command);
if ($rc != 0) {
    die("blastcompare failed. Script terminated.");
}

# Check if any alignments missed by BLAST
open(BLAST_MISSED, $ARGV[2] . $blastMissed) or die("Cannot open blast missed alignment!\n");
$tempblastMissed = <BLAST_MISSED>;
close(BLAST_MISSED);
if ($tempblastMissed !~ /\S/) {
  unlink($ARGV[2] . $blastMissed);
}

# Check if any alignments missed by BWT-SW
open(BWTSW_MISSED, $ARGV[2] . $bwtswMissed) or die("Cannot open bwtsw missed alignment!\n");
$tempbwtswMissed = <BWTSW_MISSED>;
close(BWTSW_MISSED);
if ($tempbwtswMissed !~ /\S/) {
  unlink($ARGV[2] . $bwtswMissed);
}

print "done.\n";

