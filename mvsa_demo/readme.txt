% This package contains a matlab implementation of the MVSA (hard constraint version)
% algorithm [1]. For a soft constraint version, see SISAL algorithm [2])
%
%
%--------------------------------------------------------------------------
%   Files included
%-------------------------------------------------------------------------
%
%  mvsa.m    -> MVSA algorithm [1]
%  vca.m     -> VCA algorithm [3]
%  estNoise  -> Noise estimation algorithm. See [4]
%  dataProj  -> Project data algorithm
%  hysime    -> Hysime algorithm [5]
%  USGS_1995_Library.mat -> USGS spectral library
%
%  DEMOS:
%
%  mvsa_demo_easy.m     -> Very easy problem
%  mvsa_demo_p10.m      -> medium/large number of endmembers
%  mvsa_demo_n10000.m   -> medium/large number of pixels
%  mvsa_demo_SNR20dB.m  -> low SNR
%  mvsa_demo_USGS.m     -> signatures from USGS library
%  mvsa_demo_hard.m     -> hard problem
%
%--------------------------------------------------------------------------
%    How to run
%-------------------------------------------------------------------------
%
%  Simply download the complete package to a directory and run the demos
%
% MVSA: Minimum volume simplex analysis
%
% [1] Jun Li and José M. Bioucas-Dias
%     "Minimum volume simplex analysis: A fast algorithm to unmix hyperspectral data"
%      in IEEE International Geoscience and Remote sensing Symposium
%      IGARSS’2008, Boston, USA,  2008.
%
% SISAL: Simplex identification via split augmented Lagrangian
%
% [2] J. Bioucas-Dias, "A variable splitting augmented Lagrangian approach
%     to linear spectral unmixing", in  First IEEE GRSS Workshop on
%     Hyperspectral Image and Signal Processing-WHISPERS'2009, Grenoble,
%     France,  2009. Available at http://arxiv.org/abs/0904.4635v
%
% VCA: Vertex component analysis
%
% [3] J. Nascimento and J. Bioucas-Dias, "Vertex component analysis",
%     IEEE Transactions on Geoscience and Remote Sensing, vol. 43, no. 4,
%     pp. 898-910, 2005.
%
% [4] J. Bioucas- Dias and J. Nascimento, "Hyperspectral subspace
%     identification", IEEE Transactions on Geoscience and Remote Sensing,
%     vol. 46, no. 8, pp. 2435-2445, 2008
%
%
% NOTE:  VCA (Vertex Component Analysis) is used to initialize MVSA. However,
%        VCA is a pure-pixel based algorithm and thus it is not suited to
%        the data sets herein considered.  Nevertheless, we plot VCA results,
%        to highlight the advantage of non-pure-pixel based algorithms over the
%        the pure-pixel based ones.
%
%
% Author: Jose M. Bioucas-Dias (December, 2009)
%
