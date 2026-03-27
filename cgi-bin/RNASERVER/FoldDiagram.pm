package RNASERVER::FoldDiagram;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(
    render_raw_fold_svg
    render_annotated_fold_svg
    validate_feature_ranges
);

sub render_raw_fold_svg {
    my (%args) = @_;

    my $java     = $args{java};
    my $jar      = $args{jar};
 #   my $foldout  = $args{foldout};   # optional fallback
    my $svg_out  = $args{svg_out};

    my $seq      = $args{sequence};
    my $dbn      = $args{structure};

    my $layout   = $args{layout} || 'naview';
    my $padding  = defined $args{padding} ? $args{padding} : 60;
    my $width    = defined $args{width}   ? $args{width}   : 2200;
    my $height   = defined $args{height}  ? $args{height}  : 1400;
    my $period   = defined $args{period}  ? $args{period}  : 1000; # effectively disables periodicity for raw diagram

#    my ($seq, $dbn) = parse_rnafold_output($foldout);

    my @cmd = build_varna_command(
        java    => $java,
        jar     => $jar,
        seq     => $seq,
        dbn     => $dbn,
        svg_out => $svg_out,
        layout  => $layout,
        period  => $period,
    );

    run_command(\@cmd, "raw VARNA SVG");

    fit_svg_viewbox_to_content(
        $svg_out,
        padding => $padding,
        width   => $width,
        height  => $height,
    );

    return $svg_out;
}

sub render_annotated_fold_svg {
    my (%args) = @_;

    my $java     = $args{java};
    my $jar      = $args{jar};
    #   my $foldout  = $args{foldout};   # optional fallback
    my $svg_out  = $args{svg_out};
    my $seq      = $args{sequence};
    my $dbn      = $args{structure};
    my $layout   = $args{layout}   || 'naview';
    my $features = $args{features} || {};
    my $styles   = $args{styles}   || {};
    my $period   = defined $args{period} ? $args{period} : 50;
    my $padding  = defined $args{padding} ? $args{padding} : 60;
    my $width    = defined $args{width}   ? $args{width}   : 2200;
    my $height   = defined $args{height}  ? $args{height}  : 1400;


    ## code for debugging very helpful should be kept
    # for my $feature (sort keys %$features) {
    #     my @flat = @{ $features->{$feature} || [] };
    #     next unless @flat;
    #     print "DEBUG feature $feature ranges: " . join(", ", @flat) . "\n";
    # }

    #my ($seq, $dbn) = parse_rnafold_output($foldout);
    validate_feature_ranges($features, length($seq));

    my @cmd = build_varna_command(
        java     => $java,
        jar      => $jar,
        seq      => $seq,
        dbn      => $dbn,
        svg_out  => $svg_out,
        layout   => $layout,
        features => $features,
        styles   => $styles,
        period   => $period,
    );

    run_command(\@cmd, "annotated VARNA SVG");

    fit_svg_viewbox_to_content(
        $svg_out,
        padding => $padding,
        width   => $width,
        height  => $height,
    );

    return $svg_out;
}

sub validate_feature_ranges {
    my ($features, $seq_len) = @_;

    for my $feature (sort keys %$features) {
        my @flat = @{ $features->{$feature} || [] };
        next unless @flat;

        die "Feature '$feature' has odd number of coordinates\n"
            if @flat % 2 != 0;

        my @clean;

        for (my $i = 0; $i < @flat; $i += 2) {
            my ($from, $to) = @flat[$i, $i + 1];

            die "Feature '$feature' has invalid range $from-$to\n"
                unless defined $from && defined $to
                    && $from =~ /^\d+$/ && $to =~ /^\d+$/;

            # normalize reversed order
            ($from, $to) = ($to, $from) if $from > $to;

            # skip fully outside sequence
            if ($to < 1 || $from > $seq_len) {
                warn "Skipping feature '$feature' range $from-$to: outside sequence length $seq_len\n";
                next;
            }

            # clip partially overlapping ranges
            $from = 1        if $from < 1;
            $to   = $seq_len if $to   > $seq_len;

            push @clean, $from, $to;
        }

        $features->{$feature} = \@clean;
    }
}

sub build_varna_command {
    my (%args) = @_;

    my $java     = $args{java};
    my $jar      = $args{jar};
    my $seq      = $args{seq};
    my $dbn      = $args{dbn};
    my $svg_out  = $args{svg_out};
    my $layout   = $args{layout} || 'naview';
    my $features = $args{features} || {};
    my $styles   = $args{styles}   || {};
    my $period   = $args{period};

    my @cmd = (
        $java,
        '-Djava.awt.headless=true',
        '-cp', $jar,
        'fr.orsay.lri.varna.applications.VARNAcmd',
        '-sequenceDBN',  $seq,
        '-structureDBN', $dbn,
        '-o',            $svg_out,
        '-algorithm',    $layout,
        '-bpStyle',      'lw',
        '-spaceBetweenBases', '0.5',
    );

    my @highlight_specs = build_highlight_regions($features, $styles);

    if (@highlight_specs) {
        push @cmd, '-highlightRegion', join(';', @highlight_specs);
    }

    if (defined $period) {
        push @cmd, '-periodNum', $period;
    }

    return @cmd;
}

sub build_highlight_regions {
    my ($features, $styles) = @_;

    # thius is the order, can be cahnged, would like to plug this in the mian script maybe later.
    my @feature_order = qw(
            utr  
            exons 
            polyA_signal 
            polyA_tail
            motif      
            mirna       
            trna      
            sm         
            riboswitch  
            transsplice 
            ire       
            rbp       
        );

    my @highlight_specs;

    for my $feature (@feature_order) {
        next unless exists $styles->{$feature};
        next unless exists $features->{$feature};

        my $style = $styles->{$feature};
        my $mode  = $style->{mode}  || 'region';
        my $color = $style->{color} || '#CCCCCC';
        my @flat  = @{ $features->{$feature} || [] };

        next unless @flat;
        next unless $mode eq 'region' || $mode eq 'both';

        for (my $i = 0; $i < @flat; $i += 2) {
            my ($from, $to) = @flat[$i, $i + 1];

            next unless defined $from && defined $to;
            next unless $from =~ /^\d+$/ && $to =~ /^\d+$/;

            ($from, $to) = ($to, $from) if $from > $to;

            push @highlight_specs,
                "$from-$to:fill=$color,outline=$color,radius=15";
        }
    }

    return @highlight_specs;
}

sub run_command {
    my ($cmd_ref, $label) = @_;

    print STDERR "\nRunning $label:\n";
    print STDERR join(" ", map { shell_quote($_) } @$cmd_ref), "\n\n";

    system(@$cmd_ref) == 0
        or die "ERROR: $label failed with exit code " . ($? >> 8) . "\n";
}

sub shell_quote {
    my ($s) = @_;
    return "''" if !defined($s) || $s eq '';
    $s =~ s/'/'"'"'/g;
    return "'$s'";
}

sub fit_svg_viewbox_to_content {
    my ($svg_file, %opt) = @_;

    my $padding = defined $opt{padding} ? $opt{padding} : 40;
    my $out_w   = defined $opt{width}   ? $opt{width}   : 2000;
    my $out_h   = defined $opt{height}  ? $opt{height}  : 2000;

    open my $in, '<', $svg_file or die "Cannot read $svg_file: $!";
    my @lines = <$in>;
    close $in;

    my $svg = join('', @lines);

    my (@xs, @ys);

    while ($svg =~ /\b(?:x|cx|x1|x2)="([-+]?\d*\.?\d+)"/g) {
        push @xs, $1;
    }
    while ($svg =~ /\b(?:y|cy|y1|y2)="([-+]?\d*\.?\d+)"/g) {
        push @ys, $1;
    }

    while ($svg =~ /\bpoints="([^"]+)"/g) {
        my $pts = $1;
        while ($pts =~ /([-+]?\d*\.?\d+),([-+]?\d*\.?\d+)/g) {
            push @xs, $1;
            push @ys, $2;
        }
    }

    die "No drawable coordinates found in SVG\n" unless @xs && @ys;

    my ($min_x, $max_x) = minmax(@xs);
    my ($min_y, $max_y) = minmax(@ys);

    $min_x -= $padding;
    $min_y -= $padding;
    my $vb_w = ($max_x - $min_x) + $padding;
    my $vb_h = ($max_y - $min_y) + $padding;

    for (@lines) {
        next unless /<svg\b/;

        if (/\bwidth="/) {
            s/\bwidth="[^"]*"/width="$out_w"/;
        } else {
            s/<svg\b/<svg width="$out_w"/;
        }

        if (/\bheight="/) {
            s/\bheight="[^"]*"/height="$out_h"/;
        } else {
            s/<svg\b/<svg height="$out_h"/;
        }

        if (/\bviewBox="/) {
            s/\bviewBox="[^"]*"/viewBox="$min_x $min_y $vb_w $vb_h"/;
        } else {
            s/<svg\b/<svg viewBox="$min_x $min_y $vb_w $vb_h"/;
        }

        last;
    }

    open my $out, '>', $svg_file or die "Cannot write $svg_file: $!";
    print $out @lines;
    close $out;
}

sub minmax {
    my @vals = @_;
    my $min = $vals[0];
    my $max = $vals[0];

    for my $v (@vals) {
        $min = $v if $v < $min;
        $max = $v if $v > $max;
    }

    return ($min, $max);
}

1;