#!/usr/bin/perl

my $FE = 6;
my $df = 7;
my $f;
my @matrix;
my $s;
my $SA = 0;
my $SE = 9;
my $ds = 12;
my $hi = 13;
my $FrameBitL = "FER";
my $FrameBit = "fer";
my @simchar = ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n");

my $snr;
my $snr_i;
my $rate;
my $rate_i;
my $fer_zone = 0;
my @filename;
my $SMOVE;

if ($SA < 0)
{
	$SMOVE = 0 - $SA;
	$SA = 0;
	$SE = $SE + $SMOVE;
}

for ($f = 0; $f <= $df*$FE; $f++)
{
	for ($s = $SA*$ds; $s <= $SE*$ds; $s++)
	{ 
		if ($s%$ds == 0)
		{
			if ($f%$df == 0) {$matrix[$f][$s] = ":";} else {$matrix[$f][$s] = " ";}
		}
		else
		{
			if ($f%$df == 0) {$matrix[$f][$s] = ".";} else {$matrix[$f][$s] = " ";}
		}
	}
}

my $sim = 0; 
for (@ARGV)
{
	if ($_ =~ /^([0-9]+)$/) # Iterationsanzahl
	{
		$hi = $1;
	}
	elsif ($_ =~ /^ber$/) # BER ausgeben
	{
		$FrameBitL = "BER";
		$FrameBit = "ber";
	}                                
	else
	{
	# Datei einlesen
	open(DATEI, "<$_");
	my @Inhalt = <DATEI>;
	close(DATEI);
	$filename[$sim] = $_;

	for (@Inhalt)
	{
		if ($_ =~ /<iter_var name=\"eb_n0\" value=\"(\-?[0-9]+\.?[0-9]*)\">/)
		{
			$snr = $1;
			#print "SNR $snr\n";
		}
		if ($_ =~ /<status_out name=\"error_rate_blocks\">/ and $FrameBit eq "fer")
		{
			$fer_zone = 1;
		}
		if ($_ =~ /<status_out name=\"error_rate_bits\">/ and $FrameBit eq "ber")
		{
			$fer_zone = 1;
		}
		if ($_ =~ /<\/status_out>/)
		{
			$fer_zone = 0;
		}
		if ($fer_zone == 1)
		{
			if ($_ =~ /<value iteration=\"$hi\">([0-9]+\.?[0-9]*e?-?[0-9]*)<\/value>/)
			{
				$rate = $1;
				$snr_i = abs(($snr+$SMOVE) * $ds);
				if ($rate != 0)
				{
					$rate_i = sprintf("%.0f", abs(log($rate)/log(10) * $df));
					#print "SIM $simchar[$sim]  SNR $snr $snr_i  $FrameBitL $rate $rate_i\n";
					print "SIM $simchar[$sim]  SNR $snr $FrameBitL $rate\n";
					$matrix[$rate_i][$snr_i] = $simchar[$sim];
				}
				else
				{
					#print "SIM $simchar[$sim]  SNR $snr $snr_i  $FrameBitL --\n";
					print "SIM $simchar[$sim]  SNR $snr $FrameBitL --\n";
				}
			}
		}
	}
	$sim++;
	}
}

for ($s = 0; $s < $sim; $s++)
{
	print "$simchar[$s]: $filename[$s]\n";
} 
for ($f = 0; $f <= $df*$FE; $f++)
{
	if ($f%$df == 0) {printf("%2d ", -$f/$df);} else {print "   ";}
	for ($s = $SA*$ds; $s <= $SE*$ds; $s++)
	{ 
		print $matrix[$f][$s];
	}
	print "\n";
}
print "   ";
for ($s = $SA*$ds; $s <= $SE*$ds; $s++)
{ 
	if ($s%$ds == 0) {printf("%d", $s/$ds-$SMOVE);}
	elsif ($s%$ds != 1 or ($s/$ds-$SMOVE>=0 and $s/$ds-$SMOVE<10)) {print " ";}
}
print "\n";
