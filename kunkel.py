def kunkel( protocol, params ): 
  # parse csv
  with open( 'example_order.csv' ) as fn:
    mutants = [ ]
    oligos = [ ]
    n = len( oligos )
    m = len( mutants ) # used a lot to assign wells 
    for mutant_line in fn:
        s = mutant_line.split( ',' )
        mutants += [ { 'ssDNA': s.pop( 0 ), 'label': s.pop( 0 ), 'oligos': s } ]
        oligos += s

  # reagents
  oligo_plate = protocol.ref( 'oligo_plate', None, '384-flat', storage='cold_20' ) 
  #pol_mix = protocol.ref( 'pol_mix', None, 'micro-2.0', discard=True ).well( 0 ).set_volume( '99:microliter' )
  kinase_mix = protocol.ref( 'kinase_mix', None, 'micro-2.0', discard=True ).well(0).set_volume( '99:microliter' )
  
  # synthesize, dilute
  oligo_roster = { }
  for i, j in enumerate( set( oligos ) ):
    protocol.oligosynthesize( [ { 'sequence': j, 'scale': '25nm', 'destination': oligo_plate.well( i ) } ] ) # assume these come as 100 uM stocks
    oligo_roster.update( { j: oligo_plate.well( i ) } ) # now we've got oligos and plenty of extra wells
    # 7 uL 100 uM oligo + 23 uL kinase reaction mix ( 3 uL PNK buffer, 1 uL 10 mM ATP, 1 uL NEB T4 PNK, 18 uL water per reaction )
    protocol.transfer( kinase_mix, oligo_plate.well( n + i ), '0.225:microliter' ) 
    protocol.transfer( oligo_plate.well( i ), oligo_plate.well( n + i ), '0.075:microliter' ) 
    # maybe do these two transfers outside loop w/ a distribute and a transfer 
  
  #, kinase 
  protocol.seal( oligo_plate )
  protocol.incubate( oligo_plate, "warm_37", "1:hour" )
  protocol.unseal( oligo_plate )
  
  # anneal mutants
  for i, mutant in enumerate( mutants ):
    for j, oligo in enumerate( mutant['oligos'] ):
      #print 'transfer 2:microliter from %s to well %d' % ( oligo_roster[ oligo ], i ) # translate into autoprotocol
      protocol.transfer( oligo_roster[ oligo ], oligo_plate.well( n + i ), '0.2:microliter' )

  ssDNA = protocol.ref( 'ssDNA', None, 'micro-1.5', storage="cold_20" ).well( 0 )
  ligase_mix = protocol.ref( 'ligase_mix', None, 'micro-2.0', discard=True ).well(0).set_volume( '99:microliter' )
  
  protocol.transfer( ssDNA, ligase_mix, '%f:microliter' % ( 0.200 * ( len(mutants) + 2 ) ) )
  protocol.distribute( ligase_mix, oligo_plate.wells_from( n, n + i ), '0.025:microliter' )
  protocol.seal( oligo_plate )
  ramp = [ { "cycles": 1, "steps": [{ "temperature": "%d:celsius" % ( 95 - i ), "duration": "1:minute", }] } for i in range( 70 )]
  protocol.thermocycle( oligo_plate, ramp )
  protocol.unseal( oligo_plate )

  # polymerize mutants
  protocol.distribute( pol_mix, oligo_plate.wells_from( n, n + n ), '2.2:microliter' )
  protocol.seal( oligo_plate, type="low-evaporation" )
  protocol.incubate( oligo_plate, "ambient", "90:minute" )
  protocol.unseal( oligo_plate )
  
  # transform mutants
  competent_cells = protocol.ref( 'competent_cells', None, 'micro-2.0', discard=True ).well( 0 ).set_volume( '2000:microliter' )
  protocol.distribute( competent_cells, oligo_plate.wells_from( n, n + n ), '25:microliter' )
  protocol.seal( oligo_plate )
  protocol.incubate( oligo_plate, 'cold_4', '20:minute' )
  protocol.unseal( oligo_plate )

  # plate mutants
  print 'need %d 6-well agar plates' % ( ( m * 3 ) / 6 )
  
  # pick colonies

  # grow

  # sequence

if __name__ == "__main__":
  from autoprotocol.harness import run
  run(kunkel, "Kunkel")

