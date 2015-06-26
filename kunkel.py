def kunkel(protocol, params):
    # parse csv
    mutants = [ ]
    oligos = [ ]
    for mutant_line in params['mutants'].split( '\n' ):
        s = mutant_line.split( ',' )
        mutants += [ { 'ssDNA': s.pop( 0 ), 'label': s.pop( 0 ), 'oligos': s } ]
        oligos += s

    # reagents
    oligo_plate = protocol.ref( 'oligo_plate', None, '96-pcr', storage='cold_20' ) # for synthesized oligos
    pol_mix = protocol.ref( 'pol_mix', None, 'micro-2.0', discard=True ).well( 0 ).set_volume( '%f:microliter' % ( 2.2 * ( len(mutants) + 2 ) ) )
    kinase_mix = protocol.ref( 'kinase_mix', None, 'micro-2.0', discard=True ).well(0).set_volume( '%f:microliter' % ( 5 * ( len(oligos) + 2 ) ) ) #6X
    ligase_mix = protocol.ref( 'ligase_mix', None, 'micro-2.0', discard=True ).well(0).set_volume( '%f:microliter' % ( 0.2 * ( len(mutants) + 2 ) ) )

    # synthesize oligos
    oligo_roster = { }
    for i, j in enumerate( set( oligos ) ):
        protocol.oligosynthesize( [ { 'sequence': j, 'scale': '25nm', 'destination': oligo_plate.well( i ) } ] )
        oligo_roster.update( { j: oligo_plate.well( i ) } )

    # dilute oligos to 100 uM
    # may need changing depending on how oligos come in
    protocol.dispense( oligo_plate, 'water', [ { "column": i, "volume": "25:microliter" } for i in range(0, 11 ) ] )
    # fills the fill plate since we can't fill rowwise

    # kinase oligos
    protocol.distribute( kinase_mix, oligo_plate.wells_from( 0, len(mutants) ), '5:microliter' ) # 6X kinase mix
    protocol.seal( oligo_plate )
    protocol.incubate( oligo_plate, "warm_37", "1:hour" )
    protocol.unseal( oligo_plate )

    # dilute oligos
    oligo_plate2 = protocol.ref( 'oligo_plate2', None, '96-pcr', storage='cold_20' )
    protocol.dispense( oligo_plate2, 'water', [ { "column": i, "volume": "155:microliter" } for i in range(0, 11 ) ] )
    protocol.transfer( oligo_plate.wells_from( 0, len(mutants) ), oligo_plate2.wells_from( 0, len(mutants) ), '1.5:microliter')

    # anneal mutants
    for i, mutant in enumerate( mutants ):
        for j, oligo in enumerate( mutant['oligos'] ):
            print 'transfer 2:microliter from %s to well %d' % ( oligo_roster[ oligo ], i ) # translate
            protocol.transfer( oligo_roster[ oligo ], oligo_plate2.well( i ), '2:microliter' )

    protocol.transfer( params['ssDNA'], ligase_mix, '%f:microliter' % ( 2.2 * ( len(mutants) + 2 ) ) )
    protocol.distribute( ligase_mix, oligo_plate2.wells_from( 0, len(mutants) ), '2.2:microliter' )

    protocol.seal( oligo_plate2 )
    ramp = [ { "cycles": 1, "steps": [{ "temperature": "%d:celsius" % ( 95 - i ), "duration": "1:minute", }] } for i in range( 70 )]
    protocol.thermocycle( oligo_plate2, ramp )
    protocol.unseal( oligo_plate2 )

    # polymerize mutants
    protocol.distribute( pol_mix, oligo_plate2.wells_from( 0, len(mutants) ), '2.2:microliter' )
    protocol.seal( oligo_plate2 )
    protocol.incubate( oligo_plate2, "ambient", "90:minute" )
    protocol.unseal( oligo_plate2 )

    # transform mutants
    competent_cells = protocol.ref( 'competent_cells', None, 'micro-2.0', discard=True ).well( 0 ).set_volume( '2000:microliter' )
    protocol.distribute( competent_cells, oligo_plate2.wells_from( 0, len(mutants) ), '25:microliter' )
    protocol.seal( oligo_plate2 )
    protocol.incubate( oligo_plate2, 'cold_4', '20:minute' )
    protocol.unseal( oligo_plate2 )

    # plate mutants

    # pick colonies

    # grow

    # sequence


if __name__ == '__main__':
    from autoprotocol.harness import run
    run(kunkel, "Kunkel")
