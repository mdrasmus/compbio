#!/usr/bin/env perl

use Text::BibTeX;
use strict;



sub processBibfile()
{
    my ($infile) = @_;
    my $bibfile = new Text::BibTeX::File($infile);
    my $firstComment = 1;
    
    # include standard macros, like months
    $bibfile->set_structure ('Bib');

    while (my $entry = new Text::BibTeX::Entry($bibfile)) {
        next unless $entry->parse_ok;

        # first comment in a file defines the section name
        if ($entry->type() eq "comment") {
            # comments
            if ($firstComment) {
                my $name = $entry->value();
                print "<section name=\"$name\">";
            }
            $firstComment = 0;

        } else {
            # if no comment before a paper then add a default section name
            if ($firstComment) {
                print "<section name=\"Untitled\">";
                $firstComment = 0;
            }
        
            # papers
            my $key = $entry->key();
            my $type = $entry->type();

            print <<EOF;
<paper>
    <key>$key</key>
    <type>$type</type>
EOF
            foreach my $field ($entry->fieldlist()) {
                my $value = $entry->get($field);
                $value =~ s/{|}//g;

                print "    <$field>$value</$field>\n";
            }
            
            
            my $bibtex = &makeBibTex($entry);
            
print <<EOF;

    <bibtex>$bibtex</bibtex>

</paper>

EOF

        }
    }
    
    print "</section>";
    
}


sub makeBibTex()
{
    my ($entry) = @_;
    
    my $type = $entry->type();
    my $key = $entry->key();
    
    my $str = "@"."$type\{$key,\n";
    
    foreach my $field ($entry->fieldlist()) {
        my $value = $entry->get($field);
        $str .= "    $field = {$value}\n";
    }
    
    $str .= "}\n";

    return $str;
}



sub main()
{
    #
    # print header
    #
    print <<EOF;
<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="reading.xsl"?>

<papers>
EOF

    foreach my $infile (@ARGV) {
        &processBibfile($infile);
    }


    #
    # print footer
    #
    print <<EOF;
</papers>
EOF


}

main();


