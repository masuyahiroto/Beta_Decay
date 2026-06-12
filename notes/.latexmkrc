#!/usr/bin/env perl

$lualatex = 'lualatex -synctex=1 -interaction=nonstopmode %O %S';
$pdf_mode = 4;
$max_repeat = 5;

$aux_dir = 'build';
$out_dir = '.';

$pvc_view_file_via_temporary = 0;
if ($^O eq 'linux') {
    $pdf_previewer = "xdg-open %S";
} elsif ($^O eq 'darwin') {
    $pdf_previewer = "open %S";
} else {
    $pdf_previewer = "start %S";
}

$clean_full_ext = "%R.synctex.gz";
