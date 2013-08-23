{
    if ($4 == "HO") {
	printf "%4s %4s %4s %4s ", $1, $2, $3, $4
	if ($1 < 5)
	    gain = 0.00730
	else
	    gain = 0.0115
	printf "%6.4f %6.4f %6.4f %6.4f ", gain, gain, gain, gain
	print $NF
	}
    else
    	print $0
}
