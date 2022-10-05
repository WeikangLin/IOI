# Python 3.7
# "python pyioi.py -h" for help
# developer: Liqiang Hou, Weikang Lin

import numpy as np
import pandas as pd
import os
import argparse
from scipy.interpolate import interp1d

class IOI:


# Below is to load the precalculated statistics for interpolation:
    def create_summary_stat(self, IOI,N_dim):
        b_summary = np.loadtxt('./Interpolate/b_summary_'+str(N_dim)+'.txt')
        f_P_NotSignificant=interp1d(b_summary[:,0],b_summary[:,1])   # Probability for beta<1 as function of IOI
        f_b_median=interp1d(b_summary[:,0],b_summary[:,2])           # Median of beta as function of IOI
        f_b_low1=interp1d(b_summary[:,0],b_summary[:,3])             # 68%-percentile lower limit of beta
        f_b_up1=interp1d(b_summary[:,0],b_summary[:,4])              # 68%-percentile upper limit of beta
        f_b_low2=interp1d(b_summary[:,0],b_summary[:,5])             # 95%-percentile lower limit of beta
        f_b_up2=interp1d(b_summary[:,0],b_summary[:,6])              # 95%-percentile upper limit of beta
        return f_P_NotSignificant(np.sqrt(2*IOI)), f_b_median(np.sqrt(2*IOI)), f_b_low1(np.sqrt(2*IOI)), \
                 f_b_up1(np.sqrt(2*IOI)), f_b_low2(np.sqrt(2*IOI)), f_b_up2(np.sqrt(2*IOI))

# Below is to define the ranking scheme
    def print_ranking(self, median, low1, up1):
        zones = [4, 5]   # This is quite abitrary, specify another if find more adequate.
        if median+up1<zones[0]:
            print("Moderate inconsistency")
        elif median+low1>zones[1]:
            print("Very strong inconsistency")
        elif median+low1<zones[0] and median+up1>zones[1]:
            print("Strong inconsistency")
        elif median+low1<zones[0] and median+up1>zones[0]:
            print("Moderate-to-strong inconsistency")
        elif median+low1>zones[0] and median+up1>zones[0]:
            print("Strong-to-very strong inconsitency")
        else:
            print("Ranking Not Defined!")

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
            # Drop zeros row/column
            df = df[(df.T != 0).any()]
            df = df.loc[:, (df != 0).any(axis=0)]
            return df

        else:
            print('This is not magstat or corr file')
            raise FileNotFoundError

    def extractParams(self, magstat, corr, parameters):
        try:
            ind = [np.where(magstat['parameter'] == para)[0][0] for para in parameters]
        except IndexError:
            print('\n[Error]Please input correct parameter!')
            exit()
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
        
        df = pd.DataFrame.from_dict(remainIOI, orient="index", columns=["Remain", 'Outlier'])
        if sort:
            df = df.sort_values(["Outlier"], ascending=False)
        return df


    def IOI_func(self,fn=None, out=None, params=None):
        
        if not fn or not params:
            print("No input path or no parameters are specified.")
            exit()

        dat = self.folderRead(fn)
        # Common parameter Check
        param_ = [set(x['parameter']) for x, _ in dat.values()]
        unique_params = set.union(*param_)
        common_params = set.intersection(*param_)
        print('%d parameters in total, %d common parameters, %d selected.' % (len(unique_params), len(common_params), len(params)))
    
        twoIOIs = self.calcTwoIOI(dat, params)
        print(twoIOIs)
        print()
        print('\nMulti-IOI\t', self.calcMultiIOI(dat, params))
        print('\n', self.calcRemainingIOI(dat, params))

        if out:
            if not os.path.exists(os.path.dirname(out)):
                os.makedirs(os.path.dirname(out))

            with open(out, 'w') as fp:
                fp.write('%d parameters in total, %d common parameters, %d selected.' % (len(unique_params), len(common_params), len(params)))
                fp.write('\n\n')
                fp.write(self.calcTwoIOI(dat, params).to_string())
                fp.write('\n')
                fp.write('\nMulti-IOI\t' + str(self.calcMultiIOI(dat, params)))
                fp.write('\n\n' + self.calcRemainingIOI(dat, params).to_string())
        return twoIOIs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', metavar='input', help='put input folder path here', required=True)
    parser.add_argument('-o', metavar='output', help='output file path here')
    parser.add_argument('-p', metavar='param', nargs='+', help='put all parameters you want here', required=True)

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
        if not os.path.exists(os.path.dirname(out)):
            os.makedirs(os.path.dirname(out))

        with open(out, 'w') as fp:
            fp.write('%d parameters in total, %d common parameters, %d selected.' % (len(unique_params), len(common_params), len(params)))
            fp.write('\n\n')
            fp.write(ioi.calcTwoIOI(dat, params).to_string())
            fp.write('\n')
            fp.write('\nMulti-IOI\t' + str(ioi.calcMultiIOI(dat, params)))
            fp.write('\n\n' + ioi.calcRemainingIOI(dat, params).to_string())

    
