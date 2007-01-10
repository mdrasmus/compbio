#!/usr/bin/perl -w
use strict;
use BerkeleyDB;

my $STREAM;

my $COLOR_HEADING = "#c0c0ff";
my $COLOR_ROW1    = "#ffffff";
my $COLOR_ROW2    = "#e0e0e0";
my $STYLE_BOTTOM_BORDER = "border-bottom: solid black 1px;";



sub write_results()
{
    my ($infile, $outfile, $nsamples, $prefix, $dbhash) = @_;
    
    set_stream($outfile);
    
    my $color = 1;
    my $c;
    
    begin_table();
        # column headings
        begin_row();
            begin_cell($COLOR_HEADING);
                begin_bold();
                    str("Rank");
                end_bold();
            end_cell();
            begin_cell($COLOR_HEADING);
                begin_bold();
                    str("Score");
                end_bold();
            end_cell();
            begin_cell($COLOR_HEADING);
                begin_bold();
                    str("Label");
                end_bold();
            end_cell();
            begin_cell($COLOR_HEADING);
                begin_bold();
                    str("Sample");
                end_bold();
            end_cell();
        end_row();
        
        # data
        my $row = 0;
        while (my $line = <$infile>) {        
            chomp $line;
            my ($label, $score) = split /\t/, $line;
            my $sample = $dbhash->{$label};

            $row++;            
            if ($color) {
                $c = $COLOR_ROW1;
            } else {
                $c = $COLOR_ROW2;
            }
            $color = !$color;
            
            begin_row();
                begin_cell($c,1,$STYLE_BOTTOM_BORDER);
                    str($row);
                end_cell();                
                begin_cell($c,1,$STYLE_BOTTOM_BORDER);
                    if (defined($score)) {
                        str_float($score);
                    } else {
                        str("--");
                    }
                end_cell();
                begin_cell($c,1,$STYLE_BOTTOM_BORDER);
                    if (defined($sample) && 
                        ($row <= $nsamples || $nsamples == -1))
                    {
                        $prefix =~ m!^.*/([^/]*)$!;
                        str("<a href='$1-$row.html'>");
                    }
                    
                    str($label);
                    
                    if (defined($sample) && 
                        ($row <= $nsamples || $nsamples == -1))
                    {
                        str("</a>");
                    }
                end_cell();
                begin_cell($c,1,$STYLE_BOTTOM_BORDER);
                    if (defined $sample) {
                        if (length($sample) < 200) {
                            str($sample);
                        } else {
                            str(substr($sample, 0, 200));
                            str("...");
                        }
                    } else {
                        str("--");
                    }
                end_cell();
            end_row();
            
            
            # write sample file
            if (defined $sample) {
                open(FILE, "> $prefix-$row.html");
                print FILE "<html><body>";
                print FILE $sample;
                print FILE "</body></html>";
                close(FILE);
            }
        }
    end_table();  
}


sub str_float()
{
    my $line = shift @_;
    printf $STREAM "%.5f", $line;
}

sub str() {
    my $line = shift @_;
    print $STREAM $line;
}

sub set_stream() {
    $STREAM = shift @_;
}

sub begin_table() {
    &str("<table border='0' cellspacing='0' cellpadding='3'>");
}

sub end_table() {
    &str("</table>");
}
sub begin_row() {
    &str("<tr>");
}
    
sub end_row() {
    &str("</tr>");
}
    
sub begin_cell() {
    my ($color, $ncols, $style) = @_;

    if (!defined($color)) { $color = "#ffffff"; }
    if (!defined($ncols)) { $ncols = 1; }
    if (!defined($style)) { $style = ""; }
    
    &str("<td bgcolor='$color' colspan='$ncols' valign='top' style='$style'>");
}
    
sub end_cell() {
    &str("</td>");
}

sub begin_bold() { &str("<b>"); }
sub end_bold() { &str("</b>"); }    
sub new_line() { &str("<br>"); }



sub main()
{
    if (scalar(@ARGV) < 1) {
        print "usage: write_list.pl <out file> '<title>' \n";
        print "       [-db <db file>] [-doc <rlabel> <doc file>] [-nsamples <# samples>]\n";
        print "stdin: two columns <label> [<score>]\n";
        exit 1;
    }

    # parse args
    my $out = shift @ARGV;       
    my $title = shift @ARGV;
    my $dbfile;
    my $docfile;
    my $rlabel_file;
    my $nsamples = -1;
    
    for (my $arg = shift @ARGV; scalar(@ARGV); $arg = shift @ARGV) {
        if ($arg eq "-db") {
            $dbfile = shift @ARGV;
        } elsif ($arg eq "-doc") {
            $rlabel_file = shift @ARGV;
            $docfile = shift @ARGV;            
        } elsif ($arg eq "-nsample") {
            $nsamples = shift @ARGV;
        }
    }
    
    # variables

    my %dbhash;
    my $doc;

    $out =~ m/^(.*)\.html$/;
    my $prefix = $1;

    # setup file streams
    my $outfile;
    my $infile;
    open($outfile, "> $out");
    open($infile, "<&STDIN");

    # load database
    if (defined $dbfile) {
        tie %dbhash, 'BerkeleyDB::Hash', -Filename => $dbfile, -Flags => DB_RDONLY;
    } else {
        undef %dbhash;
    }
    
    # open doc file
    if (defined $docfile) {
        my @labels;
        open(FILE, "< $rlabel_file") or die "cannot open $rlabel_file";
        @labels = <FILE>;
        close(FILE);
        
        open(FILE, "< $docfile") or die "cannot open $docfile";
        my $i = 0;
        while (my $line = <FILE>) {
            chomp $labels[$i];
            $dbhash{$labels[$i]} = $line;
            $i++;
        }
        close(FILE);
    }


    # write output
    print $outfile "<html><body><center><h3>$title</h3>";
    &write_results($infile, $outfile, $nsamples, $prefix, \%dbhash);
    print $outfile "</center></body></html>";
    
    return 0;
}


main();
