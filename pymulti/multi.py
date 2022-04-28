from scipy.io import FortranFile
import numpy as np
import pandas as pd
import io
import pymoog
from . import private
import os
from matplotlib import gridspec

pymulti_path = '{}/.pymulti/'.format(os.environ['HOME'])
multi_path = pymulti_path + 'files/multi/mul23/'

def read_atmos(path):
    '''
    Read MULTI-format atmosphere file. 
    '''
    with open(path, 'r') as file:
        content = file.readlines()
    content = [line[:-1] for line in content if line[0] != '*']
    model = {}
    model['ID'] = content[0].rstrip()
    model['d-scale type'] = content[1]
    model['logg'] = float(content[2])
    model['N'] = int(content[3])
    model['part1'] = pd.read_csv(io.StringIO('\n'.join(content[4:4+model['N']])), delim_whitespace=True, 
                                 names=['log(rhox)', 'T', 'Ne', 'Vz', 'Vturb'])
    model['part2'] = pd.read_csv(io.StringIO('\n'.join(content[4+model['N']:4+2*model['N']])), delim_whitespace=True, 
                                 names=['NH(1)', 'NH(2)', 'NH(3)', 'NH(4)', 'NH(5)', 'NP'])
    
    return model

def write_atmos(model, path, no_part2=False):
    '''
    Write MULTI-format atmosphere into specified file.
    '''
    header1 = model['part1'][:1].to_string(index=False, float_format=lambda x: "{:14.6e}".format(x)).split('\n')[0]
    header1 = '*' + header1[1:] + '\n'
    part1 = model['part1'].to_string(header=False, index=False, float_format=lambda x: "{:14.6e}".format(x)) + '\n'
    
    if not(no_part2):
        header2 = model['part2'][:1].to_string(index=False, float_format=lambda x: "{:12.4e}".format(x)).split('\n')[0]
        header2 = '*  HYDROGEN POPULATIONS\n*' + header2[1:] + '\n'
        part2 = model['part2'].to_string(header=False, index=False, float_format=lambda x: "{:12.4e}".format(x)) + '\n'
    
    with open(path, 'w') as file:
        file.writelines(model['ID'] + '\n')
        file.writelines(model['d-scale type'] + '\n')
        file.writelines('* LG G \n')
        file.writelines('  {}\n'.format(model['logg']))
        file.writelines('* NDEP \n')
        file.writelines('  {}\n'.format(model['N']))
        file.writelines('* \n')
        file.writelines(header1)
        file.writelines(part1)
        if not(no_part2):
            file.writelines(header2)
            file.writelines(part2)

def create_input_file(input_file, kwargs):
    '''
    Create MULTI input files.
    '''
    aval_vars = ['DIFF', 'ELIM1', 'ELIM2', 'QNORM', 'THIN', 'IATOM2', 'ICONV', 'IHSE', 'ILAMBD', 'IOPAC', 'ISTART', 'ISUM', 
                 'ITMAX', 'ITRAN', 'NMU', 
                 'IWABND', 'IWATMS', 'IWATOM', 'IWCHAN', 'IWDAMP', 'IWEMAX', 'IWEQW', 'IWEVEC', 'IWHEAD', 'IWHSE', 'IWLGMX', 
                 'IWLINE', 'IWLTE', 'IWN', 'IWWNIT', 'IWOPAC', 'IWRAD', 'IWSTRT', 'IWTAUQ', 'IWTEST', 'IWWMAT', 'IWJFIX', 
                 'IWARN', 'IOPACL', 'ISCAT', 'INCRAD', 'INGACC', 'ICRSW', 'IOSMET', 'EOSMET', 'IDL1', 'IDLNY', 'IDLCNT', 'IDLOPC']
    
    kwargs_default = {'DIFF':5, 'ELIM1':0.1, 'ELIM2':0.01, 'QNORM':15, 'THIN':0.1, 'IATOM2':0, 'ICONV':0,
                      'IHSE':0, 'ILAMBD':1, 'IOPAC':1, 'ISUM':0, 'ITMAX':50, 'ITRAN':0, 'NMU':5, 'IOPACL':0, 'ISCAT':1,
                      'IDL1':1, 'IDLNY':1, 'IDLCNT':1, 'IDLOPC':1}
        
    # combine kwargs_default and kwargs
    kwargs = {**kwargs_default, **kwargs}
    
    # Check if any key is not allowed
    for key in kwargs.keys():
        if key not in aval_vars:
            raise ValueError('Variable {} is not in variable list of MULTI.'.format(key))
            
    # Write the vairables into file
    output_list = []
    i = 1
    for key in aval_vars:
        if key in kwargs.keys():
            output_list.append('{}={},'.format(key, kwargs[key]))            
            if np.mod(i, 5) == 0:
                output_list.append('\n')
            i += 1
    with open(input_file, 'w') as file:
        file.writelines(output_list)
    pass

def mulrd(idl1_file):
    '''
    Read IDL1 output file from MULTI; same as the IDL routine MULRD.
    '''
    file = FortranFile(idl1_file, 'r')
    res = {}
    
    hh = 6.626176E-27
    cc = 2.99792458E+10
    iformat = 0
    ndep,nk,nline,nwide,nrad,nrfix,nmu,mq = file.read_record(dtype='int32')
    
    if ndep == 0:
        raise ValueError('ndep is 0; refer to this part of mulrd.pro.')
    nq = file.read_record(dtype='int32')
    qnorm = file.read_record(dtype='float32')
    abnd,awgt = file.read_record(dtype='float32')
    ev = file.read_record(dtype='float32')
    g = file.read_record(dtype='float32')
    ion = file.read_record(dtype='int32')
    hn3c2 = file.read_record(dtype='float32')

    if nrad > 0:
        ktrans = file.read_record(dtype='int32')
        jrad = file.read_record(dtype='int32')
        irad = file.read_record(dtype='int32')
        f = file.read_record(dtype='float32')
        iwide = file.read_record(dtype='int32')
        ga = file.read_record(dtype='float32')
        gw = file.read_record(dtype='float32')
        gq = file.read_record(dtype='float32')

    if iformat == 0:
        krad = file.read_record(dtype='int32')
    else:
        krad = np.zeros(max(irad), max(jrad))
        for kr in range(nrad):
            i = irad[kr] - 1
            j = jrad[kr] - 1
            krad[i, j] = kr + 1
            krad[j, i] = kr + 1

    if iformat == 0:
        z = file.read_record(dtype='float32')

    if nwide > 0:
        alfac = file.read_record(dtype='float32')

    hny4p = file.read_record(dtype='float32')

    if nrad > 0:
        alamb = file.read_record(dtype='float32')

    if nline > 0:
        a = file.read_record(dtype='float32')

    if iformat == 0:
        b = file.read_record(dtype='float32')
    else:
        b = np.zeros(max(irad), max(jrad))
        for kr in range(nline):
            i = irad[kr] - 1
            j = jrad[kr] - 1
            hn3c2_line = 2.*hh*cc/(alamb(kr)*1.e-8)**3
            b[i, j] = a[kr]/hn3c2_line
            b[j, i] = g[j]/g[i]*b[j,i]

    totn = file.read_record(dtype='float32')

    if nrad > 0 and iformat == 0:
        bp = file.read_record(dtype='float32').reshape(-1, ndep)

    nstar = file.read_record(dtype='float32')
    n = file.read_record(dtype='float32')
    n = n.reshape(ndep, nk)

    if iformat == 0:
        c = file.read_record(dtype='float32')
        c_dict = {}
        for j in range(1, nk+1):
            for i in range(1, nk+1):
                c_dict['{}-{}'.format(i, j)] = []
        j = 1
        i = 1
        for ele in c:
            c_dict['{}-{}'.format(i, j)].append(ele)
            if i == nk:
                i = 0
                j += 1
            if j == nk+1:
                j = 1
            i += 1
        for j in range(1, nk+1):
            for i in range(1, nk+1):
                c_dict['{}-{}'.format(i, j)] = np.array(c_dict['{}-{}'.format(i, j)])
        c = c_dict

    if nrfix > 0:
        jfx = file.read_record(dtype='float32')
        ifx = file.read_record(dtype='float32')
        ipho = file.read_record(dtype='float32')
        a0 = file.read_record(dtype='float32')
        trad = file.read_record(dtype='float32')
        itrad = file.read_record(dtype='float32')
                   
        
    dnyd = file.read_record(dtype='float32')

    if nline > 0 and iformat > 0:
        adamp = file.read_record(dtype='float32')

    _ = file.read_record(dtype='S2')
    _ = file.read_record(dtype='S2')
    atomid = file.read_record(dtype='S20')
    crout = file.read_record(dtype='S6')

    grav = file.read_record(dtype='float32')
    cmass = file.read_record(dtype='float32')
    temp = file.read_record(dtype='float32')
    nne = file.read_record(dtype='float32')
    vel = file.read_record(dtype='float32')
    tau = file.read_record(dtype='float32')
    lgtau = np.log10(tau)

    xnorm = file.read_record(dtype='float32')
    height = file.read_record(dtype='float32')
    temp_string = file.read_record(dtype='S145')
    vturb = file.read_record(dtype='float32')

    # common block catmo2
    if iformat == 0:
        bh = file.read_record(dtype='float32')
    nh = file.read_record(dtype='float32').reshape(ndep, -1)
    if iformat == 0:
        rho = file.read_record(dtype='float32')

    # common block csline
    if nrad > 0:
        qmax = file.read_record(dtype='float32')
        q0 = file.read_record(dtype='float32')
        ind = file.read_record(dtype='int32')
        diff = file.read_record(dtype='float32')
        q = file.read_record(dtype='float32')
        q = q.reshape(nrad, -1)
        if iformat == 0:
            wq = file.read_record(dtype='float32')

    wqmu = file.read_record(dtype='float32')    
    
    if nwide > 0:
        frq = file.read_record(dtype='float32')
    if nrad > 0:
        if iformat == 0:
            wphi = file.read_record(dtype='float32')
        sl = file.read_record(dtype='float32').reshape(-1, ndep)
    if nline > 0:
        weqlte = file.read_record(dtype='float32')
        weq = file.read_record(dtype='float32')
    if nrad > 0:
        if iformat == 0:
            rij = file.read_record(dtype='float32').reshape(-1, ndep)
            rji = file.read_record(dtype='float32').reshape(-1, ndep)
        if iformat != 2:
            flux = file.read_record(dtype='float32')
            flux = flux.reshape(nrad, -1)
        else:
            raise ValueError('This part is not done yet.')
        if iformat != 2:
            outint = file.read_record(dtype='float32')
            outint = outint.reshape(nrad, nmu, -1)
        else:
            raise ValueError('This part is not done yet.')

        if iformat == 0:
            cool = file.read_record(dtype='float32').reshape(-1, ndep)

    # common block cgausi
    xmu = file.read_record(dtype='float32')
    wmu = file.read_record(dtype='float32')

    # common block cconst
    ee,hh,cc,bk,em,uu,hce,hc2,hck,ek,pi = file.read_record(dtype='float32')
    
    for name in ['hh','cc','ndep','nk','nline','nwide','nrad','nrfix','nmu','mq',
                 'nq','qnorm','abnd','awgt','ev','g','ion','hn3c2',
                 'ktrans','jrad','irad','f','iwide','ga','gw','gq',
                 'krad','z','alfac','hny4p','alamb','a','b','totn',
                 'bp','nstar','n','c','jfx','ifx','ipho','a0','trad','itrad',
                 'dnyd','adamp','atmoid','crout','grav','cmass','temp','nne','vel','tau','lgtau',
                 'xnorm','height','vturb','bh','nh','rho','qmax','q0','ind','diff','q','wq','wqmu',
                 'frq','wphi','sl','weqlte','weq','rij','rji','flux','outint','cool','xmu','wmu']:
        try:
            res[name] = eval(name)
        except NameError:
            pass    
    
    return res

def extract_profile(res, kr, fitype='outint'):
    '''
    double.pro
    '''
    if fitype == 'flux':
        # two dimensional y array: flux values in y
        y2 = res['flux'][kr, 1:res['nq'][kr]+1]
    elif fitype == 'outint':
        # three dimensional y array: outint values in y
        y2 = res['outint'][kr, res['nmu']-1, 1:res['nq'][kr]+1]
        
    qn = res['qnorm'] * 1e5 / res['cc']
    
    if res['ind'][kr] == 2:
        # two sided y2, no need to reverse it
        yy = y2
        qx = res['q'][kr, 0:res['nq'][kr]]
    else:
        q2 = res['q'][kr, 0:res['nq'][kr]]
        qx = np.concatenate([-q2[::-1][:-1], q2])
        yy = np.concatenate([y2[::-1][:-1], y2])
        
    xx = - res['alamb'][kr] * qx * qn / (1 + qx*qn) + res['alamb'][kr]
    
    return xx[::-1], yy[::-1]

def initiate_dscale(model, path, N=-1, tau5000_top=-8):
    
    '''
    Initiate the depth scale. NDEP will be at least 100.
    '''
    
    column_1 = model['part1'].columns[0]
    
    if model['d-scale type'][0].lower() == 't':
        tau5000_top = model['part1'].loc[0, 'log_tau5000'] - 0.1
    
    with open(path, 'w') as file:
        file.writelines(model['ID'] + ' Equidistant' + '\n')
        file.writelines(model['d-scale type'] + '\n')
        if N < 0:
            if np.abs(N) < 100:
                file.writelines('{:.0f} {:.1f} \n'.format(-100, tau5000_top))
            else:
                file.writelines('{:.0f} {:.1f} \n'.format(-model['N'], tau5000_top))
        else:
            file.writelines('{:.0f} {:.1f} \n'.format(-N, tau5000_top))

        file.writelines('  {:12.6E} \n'.format(min(model['part1'][column_1])))
        print(min(model['part1'][column_1]), max(model['part1'][column_1]))
        file.writelines('  {:12.6E} \n'.format(max(model['part1'][column_1])))    
    pass

def create_IMINUS(iminus_df, path, SCALEI=1):
    
    '''
    Create the IMINUS file.
    '''
    
    NMUI = len(iminus_df.columns) - 1
    if NMUI != 1:
        raise ValueError('NMUI is not equal to 1. MULTI will not run for this.')
    NXL = len(iminus_df)
    
    IMINUS_content = []
    IMINUS_content.append('{}\n'.format(SCALEI))
    IMINUS_content.append('{}\n'.format(NMUI))
    IMINUS_content.append('{}\n'.format(NXL))
    IMINUS_content.append(iminus_df.to_string(header=False, index=False, float_format=lambda x: "{:14.6e}".format(x)) + '\n')
    
    with open(path, 'w') as file:
        file.writelines(IMINUS_content)

    pass

def create_dscale(model, path, tau5000_top=-8, initiate=True, N=-1, dscale=None):
    
    column_1 = model['part1'].columns[0]
    
    if model['d-scale type'][0].lower() == 't':
        tau5000_top = model['part1'].loc[0, 'log_tau5000'] - 0.1
    
    with open(path, 'w') as file:
        file.writelines(model['ID'] + '\n')
        file.writelines(model['d-scale type'] + '\n')
        if initiate:
            if N < 0:
                file.writelines(' {:.0f} {:.1f} \n'.format(-model['N'], tau5000_top))
            else:
                file.writelines(' {:.0f} {:.1f} \n'.format(-N, tau5000_top))

            file.writelines('  {:13.6E} \n'.format(min(model['part1'][column_1])))
            file.writelines('  {:13.6E} \n'.format(max(model['part1'][column_1])))    
        else:
            file.writelines(' {:.0f} {:.1f} \n'.format(len(dscale), tau5000_top))
            for num in dscale:
                file.writelines('  {:13.6E} \n'.format(num))
    pass

class multi():
    
    def __init__(self, teff, logg, m_h, T_min, Z_control_array, T_control_array, abun_change={}):
        self.multi_path = multi_path
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.T_min = T_min
        self.Z_control_array = Z_control_array
        self.T_control_array = T_control_array
        self.dscale = {}
        
        # Modify the abundance change 
        self.abun_change = abun_change
        pd.options.mode.chained_assignment = None
        self.asplund09 = private.asplund09[['Z', 'element', 'adopted']]
        self.asplund09.loc[:, 'element'] = self.asplund09.loc[:, 'element'].map(str.upper)
        self.asplund09.loc[2:, 'adopted'] = self.asplund09.loc[2:, 'adopted'] + self.m_h
        for ele in abun_change.keys():
            if ele.lower() in self.asplund09['element'].map(str.lower).values:
                print(ele, 'in element')
                self.asplund09.loc[self.asplund09['element'].map(str.lower) == ele.lower(), 'adopted'] = self.asplund09.loc[self.asplund09['element'].map(str.lower) == ele.lower(), 'adopted'] + abun_change[ele]
        self.asplund09 = self.asplund09[['element', 'adopted']]
        pd.options.mode.chained_assignment = 'warn'
        
    def create_model(self, smooth_window_len=13, min_lgrhox=-5, len_chromo=40, len_TR=10, min_lgNE=6, vturb=5, T_TR=1e5):
        # Read the Kurucz model
        model_kurucz = pymoog.model.interpolate_model(self.teff, self.logg, self.m_h, to_path=False, kurucz_format=True)
        model_kurucz = private.pd.DataFrame(model_kurucz[1], columns=['RHOX','T','P','XNE','ABROSS','ACCRAD','VTURB', 'FLXCNV','VCONV'])

        rhox_pointer = private.np.argmin(private.np.abs(model_kurucz['T'] - self.T_min))
        model_photo = model_kurucz[rhox_pointer::2]
        model_photo.reset_index(drop=True, inplace=True)
        rhox_photo = model_photo['RHOX'].values
        log_rhox_photo = private.np.log10(rhox_photo)
        TE_photo = model_photo['T'].values
        P_photo = model_photo['P'].values
        NE_photo = model_photo['XNE'].values
        
        if len(self.T_control_array) == 0:
            return rhox_photo, TE_photo, NE_photo

        # rhox
        fit_rhox_chromo = private.np.polyfit([0,len_chromo-1], [min_lgrhox, private.np.log10(model_kurucz.loc[rhox_pointer-1,'RHOX'])], 1)
        log_rhox_chromo = private.np.polyval(fit_rhox_chromo, private.np.arange(len_chromo))
        if len_TR > 0:
            fit_TR = private.np.polyfit([0, len_TR], [min_lgrhox - len_TR*1e-4, min_lgrhox], 1)
            log_rhox_TR = private.np.polyval(fit_TR, private.np.arange(len_TR))
        else:
            log_rhox_TR = private.np.array([])
        log_rhox = private.np.concatenate([log_rhox_TR, log_rhox_chromo, log_rhox_photo])
        rhox = 10**log_rhox

        T_control_array_use = self.T_control_array + [TE_photo[0]]
        Z_control_array_use = self.Z_control_array + [len_chromo]

        # TE
        TE_chromo = []
        for i in range(len(T_control_array_use)-1):
            fit_T_chormo_single = private.np.polyfit(Z_control_array_use[i:i+2], T_control_array_use[i:i+2], 1)
            TE_chromo += list(private.np.polyval(fit_T_chormo_single, private.np.arange(Z_control_array_use[i], Z_control_array_use[i+1])))
        if len_TR > 0:
            fit_TR = private.np.polyfit([0, len_TR], [T_TR, TE_chromo[0]], 1)
            TE_TR = private.np.polyval(fit_TR, private.np.arange(len_TR))
        else:
            TE_TR = private.np.array([])    
        TE = private.np.concatenate([TE_TR, TE_chromo, TE_photo])
        
        # NE
        NE_chromo = private.np.array([min(NE_photo)]*len_chromo)
        if len_TR > 0:
            res = private.np.polyfit([0, len_TR-1], [min_lgNE, private.np.log10(min(NE_photo))], 1)
            NE_TR = 10**private.np.polyval(res, np.arange(len_TR))
        else:
            NE_TR = private.np.array([])
        NE = private.np.concatenate([NE_TR, NE_chromo, NE_photo])
        
        if smooth_window_len != 0:
            TE_s = private.savgol_filter(TE, smooth_window_len, 3)
            indices = TE_s < TE
            TE_s[len_TR:len_TR+int((smooth_window_len-1)/2)][indices[len_TR:len_TR+int((smooth_window_len-1)/2)]] = TE[len_TR:len_TR+int((smooth_window_len-1)/2)][indices[len_TR:len_TR+int((smooth_window_len-1)/2)]] 
            TE = TE_s
            
            # rhox_s = 10**private.savgol_filter(log_rhox, smooth_window_len, 3)
            # if len_TR > 4:
            #     rhox = private.np.concatenate([rhox[:int(len_TR+1)], rhox_s[int(len_TR+1):]])
            # else:
            #     rhox = rhox_s
            NE_s = 10**private.savgol_filter(private.np.log10(NE), smooth_window_len, 3)
            NE = NE_s
        
        self.rhox = rhox
        self.T = TE
        self.Ne = NE
        
        Vz = np.zeros(len(rhox))
        Vturb = np.zeros(len(rhox))
        Vturb[:] = vturb
        part1 = pd.DataFrame({'log(rhox)':np.log10(rhox), 'T':TE, 'Ne':NE, 'Vz':Vz, 'Vturb':Vturb})

        # Create model file.
        mod = {}
        mod['ID'] = 'random'
        mod['d-scale type'] = 'MASS SCALE'
        mod['logg'] = self.logg
        mod['N'] = len(part1)

        mod['part1'] = part1
        self.atmos = mod
        
    def run_multi(self, ele, atmos_name, run_dir, renew_atmos=False, renew_dscale=True, initial_model=False, tau5000_top=-8, us='lus', NDEP=100, **kwargs):
        
        atmos_name_full = 'atmos.{}'.format(atmos_name)
        rstrt_name = 'rstrt.{}_{}'.format(ele, atmos_name)
        
        ele_name = private.re.findall(r"([a-z]+)[0-9]+", ele)[0]
        if len(ele_name) == 1:
            ele_name = ele_name.upper()
        elif len(ele_name) == 2:
            ele_name = ele_name[0].upper() + ele_name[1].lower()
        else:
            raise ValueErroFr('Element name ({}) format incorrect'.format(ele_name))
        
        # Captilize all the keys
        kwargs = {k.upper():v for k,v in kwargs.items()}
        
        # Create/renew atmos file
        if ele_name == 'H':
            if 'part2' in self.atmos.keys() and not initial_model:
                write_atmos(self.atmos, '{}/atmos.{}'.format(run_dir, atmos_name))
                # write_atmos(self.atmos, '{}/RSTRT'.format(run_dir))
            elif 'part2' in self.atmos.keys() and initial_model:
                print('Write model')
                write_atmos(self.atmos, '{}/atmos.{}'.format(run_dir, atmos_name))
            else:
                write_atmos(self.atmos, '{}/atmos.{}'.format(run_dir, atmos_name), no_part2=True)
        else:
            if initial_model:
                write_atmos(self.atmos, '{}/atmos.{}'.format(run_dir, atmos_name))
            else:
                write_atmos(self.atmos, '{}/RSTRT'.format(run_dir))
            
        # Create input file
        create_input_file('{}/input.{}'.format(run_dir, ele), kwargs=kwargs)
        
        # Create abund file; only those appeared in ABSDAT can be included.
        with open(self.multi_path + '/input/absdat') as file:
            absdat_content = file.readlines()
        allowed_ele_list = absdat_content[1][:-1].split()
        
        indices = []
        for allowed_ele in allowed_ele_list:
            if len(self.asplund09[self.asplund09['element'].map(str.lower) == allowed_ele.lower()].index) > 0:
                indices.append(self.asplund09[self.asplund09['element'].map(str.lower) == allowed_ele.lower()].index[0])
        self.asplund09 = self.asplund09.loc[indices].sort_index()
        
        multi_abund = []
        for i in self.asplund09.index:
            multi_abund.append('{:3} {:5.2f}\n'.format(self.asplund09.loc[i, 'element'], self.asplund09.loc[i, 'adopted']))
        with open('{}/abund'.format(run_dir), 'w') as file:
            file.writelines(multi_abund)
            
        # Create atom file
        atom_file_path = '{}/input/atom.{}'.format(self.multi_path, ele)
        print(atom_file_path)
        if os.path.isfile(atom_file_path):
            with open(atom_file_path) as file:
                atom_content = file.readlines()
            # Find the first and second line
            j = 0
            fs_line = []
            for i in range(len(atom_content)):
                if atom_content[i][0] != '*':
                    j = int(j + 1)
                else:
                    j += 1e-5
                    
                if j in [1, 2]:
                    fs_line.append(i)
                
            for ele_change in self.abun_change.keys():
                if ele_change.lower() == atom_content[fs_line[0]].split()[0].lower():
                    abun_array = np.array(atom_content[fs_line[1]].split(), dtype=float)
                    abun_array[0] += self.abun_change[ele_change]
                    atom_content[fs_line[1]] = '   {:.2f}    {:.2f}\n'.format(*list(abun_array))
            # Write to current file
            with open('{}/atom.{}'.format(run_dir, ele), 'w') as file:
                file.writelines(atom_content)
        else:
            raise ValueError('atom file not exist.')
        
        if private.re.findall(r"([a-z]+)[0-9]+", ele)[0].upper() != 'H':
            with open('{}/atom.{}'.format(run_dir, ele), 'r') as file:
                content = file.readlines()
                
            j = 0
            for i in range(len(content)):
                if content[i][0] != '*':
                    j += 1
                if j == 2:
                    break
            abund_value, awgt = np.array(private.re.findall('[-+]?[0-9]*\.?[0-9]+', content[i]), dtype=float)
            
            ele_num = private.element(ele_name).atomic_number
            if ele_num in self.abun_change.keys() and ele_name == 'He':
                abund_value = abund_value + self.abun_change[ele_num]
            elif not(ele_num in self.abun_change.keys()) and ele_name == 'He':
                pass
            else:
                if ele_num in self.abun_change.keys():
                    abund_value = abund_value + self.m_h + self.abun_change[ele_num]
                else:    
                    abund_value = abund_value + self.m_h
            content[i] = '  {:.2f}   {:.2f}\n'.format(abund_value, awgt)
            with open('{}/atom.{}'.format(run_dir, ele), 'w') as file:
                file.writelines(content)    
            
        # Judge if element is changed
        # if hasattr(self, 'ele'):
        #     if self.ele != ele and hasattr(self, 'dscale'):
        #         del self.dscale
        
        # DSCALE
        if ele_name == 'H':
            if ele_name in self.dscale.keys() and renew_dscale:
                print('Renewing dscale')
                # private.subprocess.run(['cp', run_dir+'/DSCAL2', '{}/dscale.{}_{}'.format(run_dir, ele, atmos_name)])
                with open('{}/dscale.{}_{}'.format(run_dir, ele, atmos_name), 'w') as file:
                    file.writelines(self.dscale[ele_name])
            else:
                # Initialize dscale if dscale attribute not exist.
                print('Initiating H dscale')
                initiate_dscale(self.atmos, '{}/dscale.{}_{}'.format(run_dir, ele, atmos_name), tau5000_top=tau5000_top, N=NDEP)
                with open('{}/dscale.{}_{}'.format(run_dir, ele, atmos_name), 'r') as file:
                    self.dscale[ele_name] = file.readlines()
        else:
            if ele_name in self.dscale.keys() and renew_dscale:
                print('Renewing dscale')
                # private.subprocess.run(['cp', run_dir+'/DSCAL2', '{}/dscale.{}_{}'.format(run_dir, ele, atmos_name)])
                with open('{}/dscale.{}_{}'.format(run_dir, ele, atmos_name), 'w') as file:
                    file.writelines(self.dscale[ele_name])
            else:
                # Initialize dscale if dscale attribute not exist.
                print('Initiating {} dscale using H dscale'.format(ele_name))
                with open('{}/dscale.{}_{}'.format(run_dir, ele, atmos_name), 'w') as file:
                    file.writelines(self.dscale['H'])
        
        # IOPAC: Non-LTE background opacity treatment
        if 'IOPAC' in kwargs.keys():
            if kwargs['IOPAC'] in [2,3]:
                # The run require BMET file, will move the BMET2 file in the working directory to BMET.
                os.rename('{}/BMET2'.format(run_dir, ele), '{}/BMET'.format(run_dir, ele))
        
        # Run MULTI
        private.subprocess.run(['cp', self.multi_path + '/run.sh', run_dir])
        private.subprocess.run(['cp', self.multi_path + '/input/{}'.format(atmos_name_full), run_dir])
        multi_run = private.subprocess.run(['./run.sh', ele, atmos_name, '23{}'.format(us)], stdout=private.subprocess.PIPE, cwd=run_dir)
        self.stdout = str(multi_run.stdout, encoding = "utf-8").split('\n')
        
        # Check if the run is finished normally
        with open('{}/JOBLOG'.format(run_dir), 'r') as file:
            content = file.readlines()
        if len(content) == 0:
            raise ValueError('MULTI not ended normally.')
        else:
            print('Run status: {}'.format(content[-1]))
        
        if renew_atmos:
            print('Renewing atmos')
            self.atmos = read_atmos('{}/{}'.format(run_dir, rstrt_name))
            
        if renew_dscale:
            with open(run_dir+'/DSCAL2', 'r') as file:
                self.dscale[ele_name] = file.readlines()
            
        self.ele = ele
        self.run_dir = run_dir
        
def trapez(x, y):
    
    integrand = (y + np.concatenate([y[1:], y[0:1]])) * (np.concatenate([x[1:], x[0:1]]) - x)
    return np.sum(integrand[0:len(x)-1]) / 2

def cal_center(lgtau, contri_func):
    '''
    
    '''
    diff = []
    del_lgtau = private.calculate_delta(lgtau)
    cum_sum = np.cumsum(contri_func*del_lgtau)
    diff = np.abs(np.abs(cum_sum / cum_sum[-1] - 0.5))
    return lgtau[diff.argmin()]

def cntrb(idlcnt_file, lgtau):
    '''
    Perform cntrb program in IDL: read the contribution function.
    '''
    file = FortranFile(idlcnt_file, 'r')
    res = {}
    
    ndep, nline, nrad, mq = 0, 0, 0, 0
    ndep,nline,nrad,mq = file.read_record(dtype='int32')
    
    nq = file.read_record(dtype='int32')
    cntrbi = np.zeros([nrad, mq, ndep])
    cntrbf = np.zeros([nrad, mq, ndep])
    if nline >= 0:
        cntrbr = np.zeros([nline, mq, ndep])
    xmeani = np.zeros([nrad, mq])
    xmeanf = np.zeros([nrad, mq])
    if nline >= 0:
        xmeanr = np.zeros([nline, mq])
    dummy = np.zeros(ndep)

    for kr in range(nrad):
        # print('reading contribution functions for kr={}'.format(kr))
        for ny in range(nq[kr]):
            dummy = file.read_record(dtype='float32')
            cntrbi[kr, ny, :] = dummy
            # xmeani[kr, ny] = trapez(lgtau, lgtau*dummy) / trapez(lgtau, dummy)
            xmeani[kr, ny] = cal_center(lgtau, dummy)
            dummy = file.read_record(dtype='float32')
            cntrbf[kr, ny, :] = dummy
            # xmeanf[kr, ny] = trapez(lgtau, lgtau*dummy) / trapez(lgtau, dummy)
            xmeanf[kr, ny] = cal_center(lgtau, dummy)

        if kr <= nline-1:
            for ny in range(nq[kr]):
                dummy = file.read_record(dtype='float32')
                cntrbr[kr, ny, :] = dummy
                # xmeanr[kr, ny] = trapez(lgtau, lgtau*dummy) / trapez(lgtau, dummy)
                xmeanr[kr, ny] = cal_center(lgtau, dummy)
                
        for name in ['ndep', 'nline', 'nrad', 'mq', 'nq', 'cntrbi', 'cntrbr', 'cntrbf', 'xmeani', 'xmeanr', 'xmeanf']:
            res[name] = eval(name)
            
    return res

def plotcntrb(mulrd_res, cntrb_res, kr, cntr_type='r', cntrb_range=[0, 1], t_range=[3000, 10000]):
    '''
    Plot the contribution. It have to be run after cntrb.
    '''
    wav, flux = extract_profile(mulrd_res, kr, fitype='flux')
    flux = flux / np.max(flux)
    ny = cntrb_res['nq'][kr]
    cntrb_2b = np.concatenate([cntrb_res['cntrb'+cntr_type][kr, 1:ny:, :][::-1], cntrb_res['cntrb'+cntr_type][kr, :ny, :]])
    if cntr_type != 'r':
        cntrb_2b = cntrb_2b / np.max(cntrb_2b)
    grid = np.meshgrid(mulrd_res['lgtau'], wav)

    fig = private.plt.figure(figsize=(15,4), dpi=200)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3], wspace=0.05) 
    ax0 = private.plt.subplot(gs[0])
    ax1 = private.plt.subplot(gs[1])
    ax2 = ax1.twinx()

    ax0.plot(flux, wav)
    ax0.invert_xaxis()
    
    ct = ax1.pcolor(*grid, cntrb_2b, cmap='Spectral_r', shading='auto', vmin=cntrb_range[0], vmax=cntrb_range[1], zorder=0)
    cbar = fig.colorbar(ct, pad=0.12)
    cbar.ax.set_ylabel('Contribution function')
    ax1.scatter(np.concatenate([cntrb_res['xmean'+cntr_type][kr][1:ny][::-1], 
                                cntrb_res['xmean'+cntr_type][kr][:ny]]), 
                wav, s=5, c='black', alpha=0.5, edgecolor='none')
    ax1.set_ylim(mulrd_res['alamb'][kr]-3, mulrd_res['alamb'][kr]+3)
    ax0.set_ylim(mulrd_res['alamb'][kr]-3, mulrd_res['alamb'][kr]+3)
    
    ax2.plot(mulrd_res['lgtau'], mulrd_res['temp'], c='white', zorder=100)
    ax2.set_ylim(t_range[0], t_range[1])
    ax0.set_xlabel(r'Flux')
    ax0.set_ylabel('Wavelength (A)')
    ax2.set_ylabel('T (K)')
    ax1.set_xlabel(r'$\log{\tau}$')
    ax1.get_yaxis().set_visible(False)