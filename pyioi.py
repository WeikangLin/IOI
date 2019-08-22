# Python 3.7
# 

import numpy as np
import pandas as pd
import os
import argparse

class IOI:
    _OutputFolder = 'Output'

    def readFile(self, fname):
        # Read File Data
        
        if str(fname).endswith('.margestats'):
            # .margestats file
            dat = []
            
            with open(fname, 'r') as fp:
                header = False
                while True:
                    line = fp.readline()
                    if not line:    break
                    if not header and len(line.split())>1 and line.split()[0] == 'parameter':
                        header = True
                    elif header:
                        dl = line.split()[0:3]
                        para = dl[0]
                        if para.endswith('*'):
                            para = para.rstrip('*')
                            # print('Added * paramter: %s' % para)
                        mean = float(dl[1])
                        sddev = float(dl[2])
                        dat.append([para, mean, sddev])

            df = pd.DataFrame(dat, columns=['parameter', 'mean', 'sddev'], dtype=float)
            return df

        elif str(fname).endswith('.corr'):
            # .corr File
            df = pd.read_csv(fname, delim_whitespace=True, header=None)
            return df

        else:
            print('This is not magstat or corr file')
            raise FileNotFoundError

    def extractParams(self, magstat, corr, parameters):
        ind = [np.where(magstat['parameter'] == para)[0][0] for para in parameters]
        mean, sddev = magstat.iloc[ind][['mean', 'sddev']].values.T
        cgm = np.outer(sddev, sddev)
        cov = corr.values[np.ix_(ind, ind)] * cgm
        return mean, cov
    
    
    def IOI(self, mu1, mu2, cgm1, cgm2):
        # Two experiments IOI
        return 0.5 * (mu1 - mu2) @ np.linalg.inv(cgm1 + cgm2) @ (mu1 - mu2).T


    def multiIOI(self, mus, cgms):
        # Multi experiments IOI
        n_exp = len(mus)
        n_param = len(mus[0])
        IOI_mul = 0
        mu_mul = np.zeros([n_param, ])
        F_mul = np.zeros([n_param, n_param])
        for i in range(n_exp):
            pm = np.linalg.inv(cgms[i])
            IOI_mul += mus[i] @ pm @ mus[i].T
            mu_mul += pm @ mus[i]
            F_mul += pm
        mu_mul = np.linalg.inv(F_mul) @ mu_mul
        return (IOI_mul - mu_mul.T @ F_mul @ mu_mul) / n_exp

    def folderRead(self, foldername):
        # File Check, return experiments dict
        fl = os.listdir(foldername)
        experiments = []
        for fn in fl:
            if fn.endswith('.margestats'):
                epm = fn.rstrip('.margestats')
                if epm + '.corr' not in fl:
                    print('%s not match.' % epm)
                    continue
                experiments.append(epm)
        print('%d files found, %d experiments matched.' % (len(fl)-1, len(experiments)))
        assert(len(experiments) > 1)

        dat = {}
        for exp in experiments:
            marg = self.readFile(os.path.join(foldername, exp + '.margestats'))
            corr = self.readFile(os.path.join(foldername, exp + '.corr'))
            dat[exp] = (marg, corr)
        
        return dat
    

    def calcTwoIOI(self, dat, params):
        experiments = list(dat.keys())
        n_exp = len(experiments)
        ioi = np.zeros([n_exp, n_exp])
        for i in range(n_exp):
            for j in range(i, n_exp):
                mean_i, cov_i = self.extractParams(dat[experiments[i]][0], dat[experiments[i]][1], params)
                mean_j, cov_j = self.extractParams(dat[experiments[j]][0], dat[experiments[j]][1], params)

                ioi[i, j] = self.IOI(mean_i, mean_j, cov_i, cov_j)
                ioi[j, i] = ioi[i, j]
        df = pd.DataFrame(ioi, columns=experiments, index=experiments)
        return df

    def calcMultiIOI(self, dat, params):
        experiments = list(dat.keys())
        n_exp = len(experiments)
        mus = []
        cgms = []
        for i in range(n_exp):
            mean_i, cov_i = self.extractParams(dat[experiments[i]][0], dat[experiments[i]][1], params)
            mus.append(mean_i)
            cgms.append(cov_i)
        mioi = self.multiIOI(mus, cgms)
        return mioi
        
    def calcRemainingIOI(self, dat, params, sort=True):
        experiments = list(dat.keys())
        n_exp = len(experiments)
        mioi = self.calcMultiIOI(dat, params)

        remainIOI = {}
        
        for exp in experiments:
            cdat = dat.copy()
            cdat.pop(exp)
            remain = self.calcMultiIOI(cdat, params)
            outlier = 0.5 * (n_exp * mioi - (n_exp-1) * remain)
            remainIOI[exp] = [remain, outlier]
        
        df = pd.DataFrame.from_dict(remainIOI, orient="index", columns=["Remain", 'Outlier index'])
        if sort:
            df = df.sort_values(["Outlier"], ascending=False)
        return df



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', help='put input folder path here')
    parser.add_argument('-o', help='output file path here')
    parser.add_argument('-p', nargs='+', help='put all parameters you want here')

    args = parser.parse_args()
    fn = args.r
    out = args.o
    params = args.p

    if not fn or not params:
        parser.print_help()
        exit()

    ioi = IOI()
    dat = ioi.folderRead(fn)
    # Common parameter Check
    param_ = [set(x['parameter']) for x, _ in dat.values()]
    unique_params = set.union(*param_)
    common_params = set.intersection(*param_)
    print('%d parameters in total, %d common parameters, %d selected.' % (len(unique_params), len(common_params), len(params)))
    

    print(ioi.calcTwoIOI(dat, params))
    print()
    print('\nMulti-IOI\t', ioi.calcMultiIOI(dat, params))
    print('\n', ioi.calcRemainingIOI(dat, params))

    if out:
        with open(out, 'w') as fp:
            fp.write(ioi.calcTwoIOI(dat, params).to_string())
            fp.write('\n')
            fp.write('\nMulti-IOI\t' + str(ioi.calcMultiIOI(dat, params)))
            fp.write('\n\n' + ioi.calcRemainingIOI(dat, params).to_string())

    
