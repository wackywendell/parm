Hsets='monera sharpcorrected KD aWW hydroavg'
Hscales='None max minmax zeroone'
Hjoins='arithmetic arithmeticzero geometriczero max maxzero'

#~ parallel -j32 --eta 'echo {1}-{2}-{3} ; python src/PDBtoCGYAML-hmodes.py -s {1} -r {2} -j {3} pdb/aS.pdb > startlocs/cg/aS-{1}-{2}-{3}.yml' ::: $Hsets ::: $Hscales ::: $Hjoins
#~ parallel -j32 --eta 'echo {1}-{2}-{3} ; python src/PDBtoCGYAML-sixmodes.py -s {1} -r {2} -j {3} pdb/aS.pdb > startlocs/cg/aS-six-{1}-{2}-{3}.yml' ::: $Hsets ::: $Hscales ::: $Hjoins
parallel -j4 --eta 'echo {1}-{2}-{3} ; python src/PDBtoCGYAML-hmodes.py -s {1} -r {2} -j {3} pdb/tau_abhi.pdb > startlocs/cg/tau-{1}-{2}-{3}.yml' ::: $Hsets ::: $Hscales ::: $Hjoins
parallel -j4 --eta 'echo {1}-{2}-{3} ; python src/PDBtoCGYAML-sixmodes.py -s {1} -r {2} -j {3} pdb/tau_abhi.pdb > startlocs/cg/tau-six-{1}-{2}-{3}.yml' ::: $Hsets ::: $Hscales ::: $Hjoins


#~ for Hset in $Hsets; do
#~ for Hscale in $Hscales; do
#~ for Hjoin in $Hjoins; do
    #~ f=startlocs/cg/aS-${Hset}-${Hscale}-${Hjoin}.yml
    #~ python src/PDBtoCGYAML-hmodes.py -s $Hset -r $Hscale -j $Hjoin pdb/aS.pdb > $f
    #~ f=startlocs/cg/aS-six-${Hset}-${Hscale}-${Hjoin}.yml
    #~ python src/PDBtoCGYAML-sixmodes.py -s $Hset -r $Hscale -j $Hjoin pdb/aS.pdb > $f
    #~ echo $? ${Hset}-${Hscale}-${Hjoin}
#~ done
#~ done
#~ done
