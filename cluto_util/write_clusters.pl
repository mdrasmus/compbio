#!/usr/bin/perl -w
use strict;

my $STREAM;

# external programs
my $CLUTO_STATS = "cluto_stats";
my $WRITE_LIST  = "write_list.pl";

# settings
my $COLOR_HEADING = "#c0c0ff";
my $COLOR_ROW1    = "#ffffff";
my $COLOR_ROW2    = "#e0e0e0";
my $STYLE_BOTTOM_BORDER = "border-bottom: solid black 1px;";



sub write_clustering()
{
    my ($out, $features, $sizes, $int_id, $int_per, $links) = @_;
    
    my $s = $STYLE_BOTTOM_BORDER;
    
    set_stream($out);
    
    begin_table();
        # heading
        begin_row();
            begin_cell($COLOR_HEADING);
                begin_bold();
                str("Cluster ID");
                end_bold();
            end_cell();
            begin_cell($COLOR_HEADING);
                begin_bold();            
                str("Size");
                end_bold();
            end_cell();
            begin_cell($COLOR_HEADING);
                begin_bold();            
                str("Most Descriptive Features");
                end_bold();                
            end_cell();                  
        end_row();
        
        
        # clusters
        my $color = 1;
        my $c;
        for (my $i=0; $i<scalar(@$sizes); $i++) {
            if ($color) {
                $c = $COLOR_ROW1;
            } else {
                $c = $COLOR_ROW2;
            }
            $color = !$color;
        
            begin_row();
                begin_cell($c, 1, $s);
                    if (defined($links)) {
                        str("<a href='$links->[$i]'>");
                    }
                    str("$i");
                    if (defined($links)) {
                        str("</a>");
                    }
                end_cell();
                begin_cell($c, 1, $s);
                    str("$sizes->[$i]");
                end_cell();
                begin_cell($c, 1, $s);
                    if (scalar(@$features) > 0) {
                       my @ids = @{$int_id->[$i]};
                       str("$features->[$ids[0]]");                       
                       foreach my $id (@ids[1..$#ids]) {
                          str(", $features->[$id]");
                       }
                    } else {
                        str("--");
                    }
                end_cell();
            end_row();
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

sub read_lines()
{
    my $file = shift @_;
    my @labels;
    
    open(FILE, "< $file") or die "could not open $file";
    while (my $line = <FILE>) {
        chomp $line;
        push @labels, $line;
    }
    close(FILE);
    
    return @labels;
}

sub read_part()
{
    my $file = shift @_;
    my @part;
    my @int_z;
    my @ext_z;
    
    open(FILE, "< $file") or die "could not open $file";
    while (my $line = <FILE>) {
        chomp $line;
        my ($p, $i, $e) = split / +/, $line;
        push @part, $p;
        push @int_z, $i;
        push @ext_z, $e;
    }
    close(FILE);
    
    return (\@part, \@int_z, \@ext_z);
}


sub read_stats()
{
    my ($matfile, $partfile, $nfeatures) = @_;

    my $cluster = 0;
    my @int_ids;
    my @int_percents;

    open(PIPE, "$CLUTO_STATS '$matfile' '$partfile' $nfeatures |");    
    
    while (my $line = <PIPE>) {
        my @row;
        my @row2;
        $int_ids[$cluster] = \@row;
        $int_percents[$cluster] = \@row2;

        @{$int_ids[$cluster]} = split " ", $line;
        $line = <PIPE>;
        @{$int_percents[$cluster]} = split " ", $line;
        $line = <PIPE>;
        $line = <PIPE>;
        $line = <PIPE>;

        $cluster++;
    }
    
    close(PIPE);    

    return (\@int_ids, \@int_percents);
}


sub calc_sizes
{
    my ($part) = @_;
    my @sizes;
    
    foreach my $line (@$part) {
        if (defined($sizes[$line])) {
            $sizes[$line]++;
        } else {
            $sizes[$line] = 0;
        }
    }
    
    return @sizes;
}

sub main()
{
    if (scalar(@ARGV) < 4) {
        print "usage: write_clusters.pl <data> <object labels> <part> \n";
        print "       <out file> [-db <db file>] [-doc <doc file>] [-samples <num>]\n";
        print "       [-features <labels>]\n";
        exit 1;
    }
    
    # parse args
    my $datafile = shift @ARGV;
    my $labelfile = shift @ARGV;
    my $partfile = shift @ARGV;
    my $outfile = shift @ARGV;
    my $dbfile;
    my $docfile;
    my $nsamples = 20;
    my $nfeatures = 5;
    my @features;
    
    while (@ARGV) {
        my $opt = shift @ARGV;
        
        if ($opt eq "-db") {
            $dbfile = shift @ARGV;
        } elsif ($opt eq "-samples") {
            $nsamples = shift @ARGV;
        } elsif ($opt eq "-features") {
            my $featurefile = shift @ARGV;
            @features = &read_lines($featurefile);
        } elsif ($opt eq "-doc") {
            $docfile = shift @ARGV;
        }
    }
    
    $outfile =~ m!([^/]*).html$!;
    my $prefix = $1;
    $outfile =~ m!^(.*)/[^/]+$!;
    my $dir = "";
    if (defined($1)) {
        $dir = $1;
    }
    
    
    # setup out file stream
    my $out;
    open($out, "> $outfile");
    
    # read labels
    my @labels = &read_lines($labelfile);
    my ($part, $int_z, $ext_z) = &read_part($partfile);

    # create stats
    my ($int_id, $int_per) = &read_stats($datafile, $partfile, $nfeatures);
    my @sizes = &calc_sizes($part);
    
    # setup second tier filenames
    my @links;
    my @files;
    for (my $i=0; $i<scalar(@sizes); $i++) {
        push @links, "$prefix-cluster$i.html";
        push @files, "$dir/$prefix-cluster$i.html";
    }
    
    # write main clustering summary output
    print $out "<html><body><center>";
    &write_clustering($out, \@features, \@sizes, $int_id, $int_per, \@links);
    print $out "</center></body></html>";
    close($out);
    
    # write second and third tier files
    for (my $i=0; $i<scalar(@sizes); $i++) {
        my $title = "Cluster " . ($i + 1);
        
        # open pipe
        if (defined($dbfile)) {
            open(PIPE, "| $WRITE_LIST '$files[$i]' '$title' -db '$dbfile' -nsamples '$nsamples'");
        } elsif (defined $docfile) {
            open(PIPE, "| $WRITE_LIST '$files[$i]' '$title' -doc '$labelfile' '$docfile' -nsamples '$nsamples'");        
        } else {
            open(PIPE, "| $WRITE_LIST '$files[$i]' '$title'");
        }
        
        
        # pass info to child script
        for (my $j=0; $j<scalar(@$part); $j++) {
            if ($part->[$j] == $i) {
                if (defined($int_z->[$j])) {
                    print PIPE "$labels[$j]\t$int_z->[$j]\n";
                } else {
                    print PIPE "$labels[$j]\n";
                }
            }
        }
        
        close(PIPE);
    }
    
    return 0;
}


main();
