%%% Processing script for in vivo data. %%%

%% Load image data.

load('ReconData_CSMT.mat')

SPGR_FA2_Data = csmtData{1};
SPGR_FA4_Data = csmtData{2};
SPGR_FA6_Data = csmtData{3};
SPGR_FA8_Data = csmtData{4};
SPGR_FA10_Data = csmtData{5};
SPGR_FA12_Data = csmtData{6};
SPGR_FA14_Data = csmtData{7};
SPGR_FA16_Data = csmtData{8};
SPGR_FA18_Data = csmtData{9};
SPGR_FA20_Data = csmtData{10};

SSFP_FA2_Data = csmtData{11};
SSFP_FA6_Data = csmtData{12};
SSFP_FA14_Data = csmtData{13};
SSFP_FA22_Data = csmtData{14};
SSFP_FA30_Data = csmtData{15};
SSFP_FA38_Data = csmtData{16};
SSFP_FA46_Data = csmtData{17};
SSFP_FA54_Data = csmtData{18};
SSFP_FA62_Data = csmtData{19};
SSFP_FA70_Data = csmtData{20};

SSFP0_FA2_Data = csmtData{21};
SSFP0_FA6_Data = csmtData{22};
SSFP0_FA14_Data = csmtData{23};
SSFP0_FA22_Data = csmtData{24};
SSFP0_FA30_Data = csmtData{25};
SSFP0_FA38_Data = csmtData{26};
SSFP0_FA46_Data = csmtData{27};
SSFP0_FA54_Data = csmtData{28};
SSFP0_FA62_Data = csmtData{29};
SSFP0_FA70_Data = csmtData{30};

CSMT_CatData = cat(4,abs(SPGR_FA2_Data),abs(SPGR_FA4_Data),abs(SPGR_FA6_Data),abs(SPGR_FA8_Data),abs(SPGR_FA10_Data),abs(SPGR_FA12_Data),abs(SPGR_FA14_Data),abs(SPGR_FA16_Data),abs(SPGR_FA18_Data),abs(SPGR_FA20_Data),abs(SSFP_FA2_Data),abs(SSFP_FA6_Data),abs(SSFP_FA14_Data),abs(SSFP_FA22_Data),abs(SSFP_FA30_Data),abs(SSFP_FA38_Data),abs(SSFP_FA46_Data),abs(SSFP_FA54_Data),abs(SSFP_FA62_Data),abs(SSFP_FA70_Data),abs(SSFP0_FA2_Data),abs(SSFP0_FA6_Data),abs(SSFP0_FA14_Data),abs(SSFP0_FA22_Data),abs(SSFP0_FA30_Data),abs(SSFP0_FA38_Data),abs(SSFP0_FA46_Data),abs(SSFP0_FA54_Data),abs(SSFP0_FA62_Data),abs(SSFP0_FA70_Data));
nii = make_nii(CSMT_CatData,[1.5 1.5 1.5],[0 0 0],64);
save_nii(nii,'NIFTI_CSMTData.nii');

B1Map = b1_fine;
B1_nii = make_nii(B1Map,[1.5 1.5 1.5],[0 0 0],64);
save_nii(B1_nii,'NIFTI_B1_CSMT.nii');

Data = load_nii('REGNIFTI_CSMTData.nii.gz'); Images = Data.img;

SPGR_FA2_Data = Images(:,:,:,1);
SPGR_FA4_Data = Images(:,:,:,2);
SPGR_FA6_Data = Images(:,:,:,3);
SPGR_FA8_Data = Images(:,:,:,4);
SPGR_FA10_Data = Images(:,:,:,5);
SPGR_FA12_Data = Images(:,:,:,6);
SPGR_FA14_Data = Images(:,:,:,7);
SPGR_FA16_Data = Images(:,:,:,8);
SPGR_FA18_Data = Images(:,:,:,9);
SPGR_FA20_Data = Images(:,:,:,10);

SSFP_FA2_Data = Images(:,:,:,11);
SSFP_FA6_Data = Images(:,:,:,12);
SSFP_FA14_Data = Images(:,:,:,13);
SSFP_FA22_Data = Images(:,:,:,14);
SSFP_FA30_Data = Images(:,:,:,15);
SSFP_FA38_Data = Images(:,:,:,16);
SSFP_FA46_Data = Images(:,:,:,17);
SSFP_FA54_Data = Images(:,:,:,18);
SSFP_FA62_Data = Images(:,:,:,19);
SSFP_FA70_Data = Images(:,:,:,20);

SSFP0_FA2_Data = Images(:,:,:,21);
SSFP0_FA6_Data = Images(:,:,:,22);
SSFP0_FA14_Data = Images(:,:,:,23);
SSFP0_FA22_Data = Images(:,:,:,24);
SSFP0_FA30_Data = Images(:,:,:,25);
SSFP0_FA38_Data = Images(:,:,:,26);
SSFP0_FA46_Data = Images(:,:,:,27);
SSFP0_FA54_Data = Images(:,:,:,28);
SSFP0_FA62_Data = Images(:,:,:,29);
SSFP0_FA70_Data = Images(:,:,:,30);

%% Save as NIFTI.

SPGR_FA2 = make_nii(abs(SPGR_FA2_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SPGR_FA2,'NIFTI_SPGRFA2.nii');
SPGR_FA4 = make_nii(abs(SPGR_FA4_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SPGR_FA4,'NIFTI_SPGRFA4.nii');
SPGR_FA6 = make_nii(abs(SPGR_FA6_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SPGR_FA6,'NIFTI_SPGRFA6.nii');
SPGR_FA8 = make_nii(abs(SPGR_FA8_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SPGR_FA8,'NIFTI_SPGRFA8.nii');
SPGR_FA10 = make_nii(abs(SPGR_FA10_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SPGR_FA10,'NIFTI_SPGRFA10.nii');
SPGR_FA12 = make_nii(abs(SPGR_FA12_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SPGR_FA12,'NIFTI_SPGRFA12.nii');
SPGR_FA14 = make_nii(abs(SPGR_FA14_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SPGR_FA14,'NIFTI_SPGRFA14.nii');
SPGR_FA16 = make_nii(abs(SPGR_FA16_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SPGR_FA16,'NIFTI_SPGRFA16.nii');
SPGR_FA18 = make_nii(abs(SPGR_FA18_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SPGR_FA18,'NIFTI_SPGRFA18.nii');
SPGR_FA20 = make_nii(abs(SPGR_FA20_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SPGR_FA20,'NIFTI_SPGRFA20.nii');

SSFP_FA2 = make_nii(abs(SSFP_FA2_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP_FA2,'NIFTI_SSFPFA2.nii');
SSFP_FA6 = make_nii(abs(SSFP_FA6_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP_FA6,'NIFTI_SSFPFA6.nii');
SSFP_FA14 = make_nii(abs(SSFP_FA14_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP_FA14,'NIFTI_SSFPFA14.nii');
SSFP_FA22 = make_nii(abs(SSFP_FA22_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP_FA22,'NIFTI_SSFPFA22.nii');
SSFP_FA30 = make_nii(abs(SSFP_FA30_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP_FA30,'NIFTI_SSFPFA30.nii');
SSFP_FA38 = make_nii(abs(SSFP_FA38_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP_FA38,'NIFTI_SSFPFA38.nii');
SSFP_FA46 = make_nii(abs(SSFP_FA46_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP_FA46,'NIFTI_SSFPFA46.nii');
SSFP_FA54 = make_nii(abs(SSFP_FA54_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP_FA54,'NIFTI_SSFPFA54.nii');
SSFP_FA62 = make_nii(abs(SSFP_FA62_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP_FA62,'NIFTI_SSFPFA62.nii');
SSFP_FA70 = make_nii(abs(SSFP_FA70_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP_FA70,'NIFTI_SSFPFA70.nii');

SSFP0_FA2 = make_nii(abs(SSFP0_FA2_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP0_FA2,'NIFTI_0SSFPFA2.nii');
SSFP0_FA6 = make_nii(abs(SSFP0_FA6_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP0_FA6,'NIFTI_0SSFPFA6.nii');
SSFP0_FA14 = make_nii(abs(SSFP0_FA14_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP0_FA14,'NIFTI_0SSFPFA14.nii');
SSFP0_FA22 = make_nii(abs(SSFP0_FA22_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP0_FA22,'NIFTI_0SSFPFA22.nii');
SSFP0_FA30 = make_nii(abs(SSFP0_FA30_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP0_FA30,'NIFTI_0SSFPFA30.nii');
SSFP0_FA38 = make_nii(abs(SSFP0_FA38_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP0_FA38,'NIFTI_0SSFPFA38.nii');
SSFP0_FA46 = make_nii(abs(SSFP0_FA46_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP0_FA46,'NIFTI_0SSFPFA46.nii');
SSFP0_FA54 = make_nii(abs(SSFP0_FA54_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP0_FA54,'NIFTI_0SSFPFA54.nii');
SSFP0_FA62 = make_nii(abs(SSFP0_FA62_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP0_FA62,'NIFTI_0SSFPFA62.nii');
SSFP0_FA70 = make_nii(abs(SSFP0_FA70_Data),[1.5 1.5 1.5],[0 0 0],64);
save_nii(SSFP0_FA70,'NIFTI_0SSFPFA70.nii');

%% Import registered NIFTIs.

SPGRFA2 = load_nii('NIFTI_SPGRFA2.nii'); SPGRFA2_Image = SPGRFA2.img;
SPGRFA4 = load_nii('NIFTI_SPGRFA4.nii'); SPGRFA4_Image = SPGRFA4.img;
SPGRFA6 = load_nii('NIFTI_SPGRFA6.nii'); SPGRFA6_Image = SPGRFA6.img;
SPGRFA8 = load_nii('NIFTI_SPGRFA8.nii'); SPGRFA8_Image = SPGRFA8.img;
SPGRFA10 = load_nii('NIFTI_SPGRFA10.nii'); SPGRFA10_Image = SPGRFA10.img;
SPGRFA12 = load_nii('NIFTI_SPGRFA12.nii'); SPGRFA12_Image = SPGRFA12.img;
SPGRFA14 = load_nii('NIFTI_SPGRFA14.nii'); SPGRFA14_Image = SPGRFA14.img;
SPGRFA16 = load_nii('NIFTI_SPGRFA16.nii'); SPGRFA16_Image = SPGRFA16.img;
SPGRFA18 = load_nii('NIFTI_SPGRFA18.nii'); SPGRFA18_Image = SPGRFA18.img;
SPGRFA20 = load_nii('NIFTI_SPGRFA20.nii'); SPGRFA20_Image = SPGRFA20.img;

SSFPFA2 = load_nii('NIFTI_SSFPFA2.nii'); SSFPFA2_Image = SSFPFA2.img;
SSFPFA6 = load_nii('NIFTI_SSFPFA6.nii'); SSFPFA6_Image = SSFPFA6.img;
SSFPFA14 = load_nii('NIFTI_SSFPFA14.nii'); SSFPFA14_Image = SSFPFA14.img;
SSFPFA22 = load_nii('NIFTI_SSFPFA22.nii'); SSFPFA22_Image = SSFPFA22.img;
SSFPFA30 = load_nii('NIFTI_SSFPFA30.nii'); SSFPFA30_Image = SSFPFA30.img;
SSFPFA38 = load_nii('NIFTI_SSFPFA38.nii'); SSFPFA38_Image = SSFPFA38.img;
SSFPFA46 = load_nii('NIFTI_SSFPFA46.nii'); SSFPFA46_Image = SSFPFA46.img;
SSFPFA54 = load_nii('NIFTI_SSFPFA54.nii'); SSFPFA54_Image = SSFPFA54.img;
SSFPFA62 = load_nii('NIFTI_SSFPFA62.nii'); SSFPFA62_Image = SSFPFA62.img;
SSFPFA70 = load_nii('NIFTI_SSFPFA70.nii'); SSFPFA70_Image = SSFPFA70.img;

SSFPFA2PC0 = load_nii('NIFTI_0SSFPFA2.nii'); SSFPFA2PC0_Image = SSFPFA2PC0.img;
SSFPFA6PC0 = load_nii('NIFTI_0SSFPFA6.nii'); SSFPFA6PC0_Image = SSFPFA6PC0.img;
SSFPFA14PC0 = load_nii('NIFTI_0SSFPFA14.nii'); SSFPFA14PC0_Image = SSFPFA14PC0.img;
SSFPFA22PC0 = load_nii('NIFTI_0SSFPFA22.nii'); SSFPFA22PC0_Image = SSFPFA22PC0.img;
SSFPFA30PC0 = load_nii('NIFTI_0SSFPFA30.nii'); SSFPFA30PC0_Image = SSFPFA30PC0.img;
SSFPFA38PC0 = load_nii('NIFTI_0SSFPFA38.nii'); SSFPFA38PC0_Image = SSFPFA38PC0.img;
SSFPFA46PC0 = load_nii('NIFTI_0SSFPFA46.nii'); SSFPFA46PC0_Image = SSFPFA46PC0.img;
SSFPFA54PC0 = load_nii('NIFTI_0SSFPFA54.nii'); SSFPFA54PC0_Image = SSFPFA54PC0.img;
SSFPFA62PC0 = load_nii('NIFTI_0SSFPFA62.nii'); SSFPFA62PC0_Image = SSFPFA62PC0.img;
SSFPFA70PC0 = load_nii('NIFTI_0SSFPFA70.nii'); SSFPFA70PC0_Image = SSFPFA70PC0.img;

B1_Map = load_nii('REGNIFTI_B1_CSMT.nii.gz'); B1_Image = B1_Map.img;
%B1_Map = load_nii('NIFTIREG_B1.nii.gz'); B1_Image = B1_Map.img;

WM_Mask = load_nii('BET2_NIFTI_SPGRFA20_pve_2.nii'); WM_Mask_Image = WM_Mask.img;

%%% These files are then passed into InVivo_Complete.m. %%%
