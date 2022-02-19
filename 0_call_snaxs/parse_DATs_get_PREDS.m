function parse_DATs_get_PREDS()


folder_DATs = '/home/dan/Research/data/SNS/ARCS/BKBO/data_COND/dat_T_300K/'
folder_preds = './preds/'
PAR = auto_PAR(EXP_BKBO_ARCS);

Ei = PAR.EXP.efixed
chopfreq = PAR.EXP.chopfreq

filelist_DATs = {
            'condensed_EXTRACTED__HKL__H_0.000__K_0.000__L_0.000__T_300K____20140525T153052_.mat',
            'condensed_EXTRACTED__HKL__H_0.100__K_0.000__L_0.000__T_300K____20140525T172735_.mat',
            'condensed_EXTRACTED__HKL__H_0.100__K_0.100__L_0.000__T_300K____20140525T162839_.mat',
            'condensed_EXTRACTED__HKL__H_0.100__K_0.100__L_0.100__T_300K____20140525T154109_.mat',
            'condensed_EXTRACTED__HKL__H_0.200__K_0.000__L_0.000__T_300K____20140525T173631_.mat',
            'condensed_EXTRACTED__HKL__H_0.200__K_0.100__L_0.000__T_300K____20140525T182935_.mat',
            'condensed_EXTRACTED__HKL__H_0.200__K_0.100__L_0.100__T_300K____20140525T223257_.mat',
            'condensed_EXTRACTED__HKL__H_0.200__K_0.200__L_0.000__T_300K____20140525T164324_.mat',
            'condensed_EXTRACTED__HKL__H_0.200__K_0.200__L_0.100__T_300K____20140526T001221_.mat',
            'condensed_EXTRACTED__HKL__H_0.200__K_0.200__L_0.200__T_300K____20140525T155140_.mat',
            'condensed_EXTRACTED__HKL__H_0.300__K_0.000__L_0.000__T_300K____20140525T174525_.mat',
            'condensed_EXTRACTED__HKL__H_0.300__K_0.100__L_0.000__T_300K____20140525T185750_.mat',
            'condensed_EXTRACTED__HKL__H_0.300__K_0.100__L_0.100__T_300K____20140525T230043_.mat',
            'condensed_EXTRACTED__HKL__H_0.300__K_0.200__L_0.000__T_300K____20140525T201325_.mat',
            'condensed_EXTRACTED__HKL__H_0.300__K_0.200__L_0.100__T_300K____20140526T010457_.mat',
            'condensed_EXTRACTED__HKL__H_0.300__K_0.200__L_0.200__T_300K____20140526T055504_.mat',
            'condensed_EXTRACTED__HKL__H_0.300__K_0.300__L_0.000__T_300K____20140525T165740_.mat',
            'condensed_EXTRACTED__HKL__H_0.300__K_0.300__L_0.100__T_300K____20140526T025622_.mat',
            'condensed_EXTRACTED__HKL__H_0.300__K_0.300__L_0.200__T_300K____20140526T070517_.mat',
            'condensed_EXTRACTED__HKL__H_0.300__K_0.300__L_0.300__T_300K____20140525T160144_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.000__L_0.000__T_300K____20140525T175417_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.100__L_0.000__T_300K____20140525T192650_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.100__L_0.100__T_300K____20140525T232759_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.200__L_0.000__T_300K____20140525T204205_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.200__L_0.100__T_300K____20140526T015731_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.200__L_0.200__T_300K____20140526T062206_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.300__L_0.000__T_300K____20140525T212914_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.300__L_0.100__T_300K____20140526T034911_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.300__L_0.200__T_300K____20140526T075615_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.300__L_0.300__T_300K____20140526T101855_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.400__L_0.000__T_300K____20140525T171221_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.400__L_0.100__T_300K____20140526T044656_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.400__L_0.200__T_300K____20140526T090015_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.400__L_0.300__T_300K____20140526T110149_.mat',
            'condensed_EXTRACTED__HKL__H_0.400__K_0.400__L_0.400__T_300K____20140525T161130_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.000__L_0.000__T_300K____20140525T180025_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.100__L_0.000__T_300K____20140525T194433_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.100__L_0.100__T_300K____20140525T234430_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.200__L_0.000__T_300K____20140525T210015_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.200__L_0.100__T_300K____20140526T022836_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.200__L_0.200__T_300K____20140526T063817_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.300__L_0.000__T_300K____20140525T214741_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.300__L_0.100__T_300K____20140526T042012_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.300__L_0.200__T_300K____20140526T082652_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.300__L_0.300__T_300K____20140526T103502_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.400__L_0.000__T_300K____20140525T220531_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.400__L_0.100__T_300K____20140526T051725_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.400__L_0.200__T_300K____20140526T093453_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.400__L_0.300__T_300K____20140526T113252_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.400__L_0.400__T_300K____20140526T115932_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.500__L_0.000__T_300K____20140525T171839_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.500__L_0.100__T_300K____20140526T052753_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.500__L_0.200__T_300K____20140526T094745_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.500__L_0.300__T_300K____20140526T114331_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.500__L_0.400__T_300K____20140526T120955_.mat',
            'condensed_EXTRACTED__HKL__H_0.500__K_0.500__L_0.500__T_300K____20140525T161403_.mat',
        };


pred_files = cell();
for ifile = 1:length(filelist_DATs)

    %% load DAT
    fname = filelist_DATs{ifile};
    DAT = load_DAT_file([folder_DATs, fname]);
    Q_hkl = cast(DAT.Q_hkl, 'double');

    %% generate predictions
    [centers, widths, heights] = query_snaxs_startvars(PAR, Q_hkl);

    %% save for later
    HKL_start = strfind(fname,'__HKL_');
    tag_end = strfind(fname, '____');
    label = fname(HKL_start+7 : tag_end-1);
    pfile = ["PREDS___", label, '___.mat']
    pred_files{ifile} = pfile;

    save([folder_preds, pfile], 'centers', 'widths', 'heights', 'Ei', 'chopfreq');
end