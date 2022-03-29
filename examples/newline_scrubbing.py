import alhambra.mixes as am

if __name__ == '__main__':

    s1 = am.Strand('s1', '100 uM', plate='plate 1', well='A1')
    s2 = am.Strand('s2', '100 uM', plate='plate 1', well='A2')
    mix1 = am.Mix(name='my mix',
                  actions=[am.MultiFixedConcentration(components=[s1, s2], fixed_concentration='10uM')],
                  fixed_total_volume='100 uL',
                  test_tube_name='my\nmix\n1')
    # print(mix1)
    # print('*' * 99)
    mix2 = am.Mix(name='my mix',
                  actions=[am.MultiFixedConcentration(components=[mix1, s1, s2], fixed_concentration='1uM')],
                  fixed_total_volume='100 uL',
                  test_tube_name='my\nmix\n2')
    t = mix2.table()
    print(t)