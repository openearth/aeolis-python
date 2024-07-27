import os
import re
import pandas
from datetime import datetime, timedelta

def read_uurgeg(fname):

    block = None
    blocks = {}
    legend = {}
    header = []
    data = []

    if os.path.exists(fname):
        with open(fname, 'r') as fp:
            for line in fp:
                
                # information blocks
                m = re.search('^\s*$', line)
                if m is not None:
                    block = None
                    continue
                
                m = re.search('^(.+):\s*$', line)
                if m is not None:
                    block = m.groups()[0]
                    blocks[block] = []
                    continue
                    
                if block is not None:
                    blocks[block].append(re.sub('\s+$', '', line))
                    continue
        
                # legend
                m = re.search('^([A-Z0-9]+)\s+=\s+(.+?) / (.+?)\s*$', line)
                if m is not None:
                    key, nl, en = m.groups()
                    legend[key] = {'nl':nl, 'en':en}
                
                # header
                m = re.search('^#\s*([A-Z0-9,\s]+?)\s*$', line)
                if m is not None:
                    header = re.split('\s*,\s*', m.groups()[0])
                    
                # data
                m = re.search('^\s*([0-9,\s]+?)\s*$', line)
                if m is not None:
                    data.append([int(x) if len(x) > 0 else None for x in re.split('\s*,\s*', m.groups()[0])])
                    
        df = pandas.DataFrame(data, columns=header)
        idx = [datetime.strptime('%08d%02d' % p, '%Y%m%d%H') + timedelta(minutes=30) for p in zip(df['YYYYMMDD'], df['HH']-1)]
        df = pandas.DataFrame(data, columns=header, index=idx)
        df.index = df.index.tz_localize('UTC')

        return df, legend
    else:
        raise IOError('File not found: %s' % fname)
