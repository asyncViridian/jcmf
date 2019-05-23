# Arguments
# $1 = directory to list scores in

(grep "pair posterior" $1/scores/* | sort -n -k4) > $1/pairPosterior.listscores
echo "Wrote list of pair posterior scores to $1/pairPosterior.listscores"

(grep "RNA posterior" $1/scores/* | sort -n -k4) > $1/rnaPosteriorScores.listscores
echo "Wrote list of RNA posterior scores to  $1/rnaPosteriorScores.listscores"
