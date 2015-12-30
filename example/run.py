from aeolis.model import AeoLiSWrapper

for tide_var in [False, True]:
    for wind_var in [False, True]:
        for bar in [False, True]:

            fname = 'aeolis_tide%d_wind%d_bar%d' % (tide_var, wind_var, bar)
            
            model = AeoLiSWrapper('aeolis.txt')
            model.set_configfile('%s.txt' % fname)
            model.set_params(output_file='%s.nc' % fname)
            model.set_params(tide_file='tide_sine.txt' if tide_var else 'tide.txt')
            model.set_params(wind_file='wind_random3.txt' if wind_var else 'wind.txt')

            if bar:
                model.run(callback='bar.py:add_bar')
            else:
                model.run()
