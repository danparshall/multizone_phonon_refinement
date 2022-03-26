function start_here();


%folder_snaxs = "~/code/snaxs/0_subroutines/"
folder_snaxs = "/home/dan/code/snaxs/";
folder_EXP = "/home/dan/code/bkbo/m-files/"; % EXP_BKBO_ARCS

path_EXTR = 'condensed_EXTRACTED__HKL__H_0.000__K_0.000__L_0.000__T_low____20140521T213302_.mat'
path_EXTR = '/home/dan/code/bkbo/data/condensed_EXTRACTED__HKL__H_0.500__K_0.500__L_0.500__T_300K____20140525T161403_.mat'
path_EXTR = '~/code/bkbo/data/condensed_EXTRACTED__HKL__H_0.300__K_0.000__L_0.000__T_300K____20140525T174525_.mat'
path_EXTR = '/home/dan/code/bkbo/data/condensed_EXTRACTED__HKL__H_0.000__K_0.000__L_0.000__T_300K____20140525T153052_.mat'

path_EXTR = '~/code/bkbo/data/condensed_EXTRACTED__HKL__H_0.500__K_0.000__L_0.000__T_300K____20140525T180025_.mat'


addpath(folder_snaxs);
addpath([folder_snaxs, '0_subroutines/']);
%addpath(folder_EXP);
PAR = auto_PAR(EXP_BKBO_ARCS);
DAT = load_DAT_file(path_EXTR);


Q_hkl = cast(DAT.Q_hkl, 'double');		% Gamma point is all integer, can cause crash, so cast explicitly as double

if 0
	% === initialize arrays, index variables ===
	[XTAL,EXP,INFO,PLOT,DATA,VECS] = params_fetch(PAR);
	INFO.Q = Q_hkl(1, :);
	DATA = make_DATA(PAR);
	INFO.Q_npts = size(Q_hkl,1);
	[unique_tau, cellarray_qs, Q_hkl, Q_delta] = generate_tau_q_from_Q(PAR,Q_hkl);


	% === generate VECS (including structure factor) ===
	PAR = params_update(XTAL,EXP,INFO,PLOT,DATA,VECS);
	PAR = simulate_multiQ(PAR, Q_hkl);

	TOL_ENG = 0.0000001;
	assert(isclose(PAR.VECS.energies(:,1), mean(PAR.VECS.energies, 2), TOL_ENG));  % sanity-check that energies truly are the same at all Q


	% === get energy and widths ===
	[centers, iWid] = unique(PAR.VECS.energies);
	widths = PAR.VECS.phWidths(iWid);	% each center has a unique intrinsic linewidth


	% === calculate height from structure factor ===
	hts_array = calc_height_multiQ(PAR, Q_hkl);
	height_nan = ~isnan(hts_array);


	heights = zeros(length(centers), INFO.Q_npts);
	for cen = 1:length(centers)
		idxCen = find(PAR.VECS.energies(:, 1) == centers(cen));
		heights(cen, :) = sum(hts_array(idxCen, :), 1); %/sum(height_nan(idxCen, :), 1);
	end
else
	[centers, widths, heights] = query_snaxs_startvars(PAR, Q_hkl);

end

disp("")
disp(["Saving predictions ..."])
disp(['    ... ', num2str(length(Q_hkl)) ' Q-points'])
PRED.Q_hkl = Q_hkl;
disp(['    ... ', num2str(length(centers)) ' phonon centers'])
for i_cen = 1:length(centers)
    disp(['        ', num2str(centers(i_cen), '%05.3f')])
end
PRED.centers = centers;
PRED.widths = widths;
PRED.heights = heights;
save('PREDS_test.mat', 'Q_hkl', 'centers', 'widths', 'heights');