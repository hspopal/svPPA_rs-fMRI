%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/haroonpopal/Google_Drive/svPPA_rs-fMRI/scan_data/n89_subject.list/allrois.n89_subject.list.z.summary.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2019/10/11 15:16:42

%% Initialize variables.
filename = '/Users/haroonpopal/Google_Drive/svPPA_rs-fMRI/scan_data/n89_subject.list/allrois.n89_subject.list.z.summary.csv';
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
%   column31: double (%f)
%	column32: double (%f)
%   column33: double (%f)
%	column34: double (%f)
%   column35: double (%f)
%	column36: double (%f)
%   column37: double (%f)
%	column38: double (%f)
%   column39: double (%f)
%	column40: double (%f)
%   column41: double (%f)
%	column42: double (%f)
%   column43: double (%f)
%	column44: double (%f)
%   column45: double (%f)
%	column46: double (%f)
%   column47: double (%f)
%	column48: double (%f)
%   column49: double (%f)
%	column50: double (%f)
%   column51: double (%f)
%	column52: double (%f)
%   column53: double (%f)
%	column54: double (%f)
%   column55: double (%f)
%	column56: double (%f)
%   column57: double (%f)
%	column58: double (%f)
%   column59: double (%f)
%	column60: double (%f)
%   column61: double (%f)
%	column62: double (%f)
%   column63: double (%f)
%	column64: double (%f)
%   column65: double (%f)
%	column66: double (%f)
%   column67: double (%f)
%	column68: double (%f)
%   column69: double (%f)
%	column70: double (%f)
%   column71: double (%f)
%	column72: double (%f)
%   column73: double (%f)
%	column74: double (%f)
%   column75: double (%f)
%	column76: double (%f)
%   column77: double (%f)
%	column78: double (%f)
%   column79: double (%f)
%	column80: double (%f)
%   column81: double (%f)
%	column82: double (%f)
%   column83: double (%f)
%	column84: double (%f)
%   column85: double (%f)
%	column86: double (%f)
%   column87: double (%f)
%	column88: double (%f)
%   column89: double (%f)
%	column90: double (%f)
%   column91: double (%f)
%	column92: double (%f)
%   column93: double (%f)
%	column94: double (%f)
%   column95: double (%f)
%	column96: double (%f)
%   column97: double (%f)
%	column98: double (%f)
%   column99: double (%f)
%	column100: double (%f)
%   column101: double (%f)
%	column102: double (%f)
%   column103: double (%f)
%	column104: double (%f)
%   column105: double (%f)
%	column106: double (%f)
%   column107: double (%f)
%	column108: double (%f)
%   column109: double (%f)
%	column110: double (%f)
%   column111: double (%f)
%	column112: double (%f)
%   column113: double (%f)
%	column114: double (%f)
%   column115: double (%f)
%	column116: double (%f)
%   column117: double (%f)
%	column118: double (%f)
%   column119: double (%f)
%	column120: double (%f)
%   column121: double (%f)
%	column122: double (%f)
%   column123: double (%f)
%	column124: double (%f)
%   column125: double (%f)
%	column126: double (%f)
%   column127: double (%f)
%	column128: double (%f)
%   column129: double (%f)
%	column130: double (%f)
%   column131: double (%f)
%	column132: double (%f)
%   column133: double (%f)
%	column134: double (%f)
%   column135: double (%f)
%	column136: double (%f)
%   column137: double (%f)
%	column138: double (%f)
%   column139: double (%f)
%	column140: double (%f)
%   column141: double (%f)
%	column142: double (%f)
%   column143: double (%f)
%	column144: double (%f)
%   column145: double (%f)
%	column146: double (%f)
%   column147: double (%f)
%	column148: double (%f)
%   column149: double (%f)
%	column150: double (%f)
%   column151: double (%f)
%	column152: double (%f)
%   column153: double (%f)
%	column154: double (%f)
%   column155: double (%f)
%	column156: double (%f)
%   column157: double (%f)
%	column158: double (%f)
%   column159: double (%f)
%	column160: double (%f)
%   column161: double (%f)
%	column162: double (%f)
%   column163: double (%f)
%	column164: double (%f)
%   column165: double (%f)
%	column166: double (%f)
%   column167: double (%f)
%	column168: double (%f)
%   column169: double (%f)
%	column170: double (%f)
%   column171: double (%f)
%	column172: double (%f)
%   column173: double (%f)
%	column174: double (%f)
%   column175: double (%f)
%	column176: double (%f)
%   column177: double (%f)
%	column178: double (%f)
%   column179: double (%f)
%	column180: double (%f)
%   column181: double (%f)
%	column182: double (%f)
%   column183: double (%f)
%	column184: double (%f)
%   column185: double (%f)
%	column186: double (%f)
%   column187: double (%f)
%	column188: double (%f)
%   column189: double (%f)
%	column190: double (%f)
%   column191: double (%f)
%	column192: double (%f)
%   column193: double (%f)
%	column194: double (%f)
%   column195: double (%f)
%	column196: double (%f)
%   column197: double (%f)
%	column198: double (%f)
%   column199: double (%f)
%	column200: double (%f)
%   column201: double (%f)
%	column202: double (%f)
%   column203: double (%f)
%	column204: double (%f)
%   column205: double (%f)
%	column206: double (%f)
%   column207: double (%f)
%	column208: double (%f)
%   column209: double (%f)
%	column210: double (%f)
%   column211: double (%f)
%	column212: double (%f)
%   column213: double (%f)
%	column214: double (%f)
%   column215: double (%f)
%	column216: double (%f)
%   column217: double (%f)
%	column218: double (%f)
%   column219: double (%f)
%	column220: double (%f)
%   column221: double (%f)
%	column222: double (%f)
%   column223: double (%f)
%	column224: double (%f)
%   column225: double (%f)
%	column226: double (%f)
%   column227: double (%f)
%	column228: double (%f)
%   column229: double (%f)
%	column230: double (%f)
%   column231: double (%f)
%	column232: double (%f)
%   column233: double (%f)
%	column234: double (%f)
%   column235: double (%f)
%	column236: double (%f)
%   column237: double (%f)
%	column238: double (%f)
%   column239: double (%f)
%	column240: double (%f)
%   column241: double (%f)
%	column242: double (%f)
%   column243: double (%f)
%	column244: double (%f)
%   column245: double (%f)
%	column246: double (%f)
%   column247: double (%f)
%	column248: double (%f)
%   column249: double (%f)
%	column250: double (%f)
%   column251: double (%f)
%	column252: double (%f)
%   column253: double (%f)
%	column254: double (%f)
%   column255: double (%f)
%	column256: double (%f)
%   column257: double (%f)
%	column258: double (%f)
%   column259: double (%f)
%	column260: double (%f)
%   column261: double (%f)
%	column262: double (%f)
%   column263: double (%f)
%	column264: double (%f)
%   column265: double (%f)
%	column266: double (%f)
%   column267: double (%f)
%	column268: double (%f)
%   column269: double (%f)
%	column270: double (%f)
%   column271: double (%f)
%	column272: double (%f)
%   column273: double (%f)
%	column274: double (%f)
%   column275: double (%f)
%	column276: double (%f)
%   column277: double (%f)
%	column278: double (%f)
%   column279: double (%f)
%	column280: double (%f)
%   column281: double (%f)
%	column282: double (%f)
%   column283: double (%f)
%	column284: double (%f)
%   column285: double (%f)
%	column286: double (%f)
%   column287: double (%f)
%	column288: double (%f)
%   column289: double (%f)
%	column290: double (%f)
%   column291: double (%f)
%	column292: double (%f)
%   column293: double (%f)
%	column294: double (%f)
%   column295: double (%f)
%	column296: double (%f)
%   column297: double (%f)
%	column298: double (%f)
%   column299: double (%f)
%	column300: double (%f)
%   column301: double (%f)
%	column302: double (%f)
%   column303: double (%f)
%	column304: double (%f)
%   column305: double (%f)
%	column306: double (%f)
%   column307: double (%f)
%	column308: double (%f)
%   column309: double (%f)
%	column310: double (%f)
%   column311: double (%f)
%	column312: double (%f)
%   column313: double (%f)
%	column314: double (%f)
%   column315: double (%f)
%	column316: double (%f)
%   column317: double (%f)
%	column318: double (%f)
%   column319: double (%f)
%	column320: double (%f)
%   column321: double (%f)
%	column322: double (%f)
%   column323: double (%f)
%	column324: double (%f)
%   column325: double (%f)
%	column326: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
SUBJECTS = dataArray{:, 1};
rV1p_lV1p = dataArray{:, 2};
lV1c_lV1p = dataArray{:, 3};
lV1c_rV1p = dataArray{:, 4};
rV1c_lV1p = dataArray{:, 5};
rV1c_rV1p = dataArray{:, 6};
rV1c_lV1c = dataArray{:, 7};
lMT_lV1p = dataArray{:, 8};
lMT_rV1p = dataArray{:, 9};
lMT_lV1c = dataArray{:, 10};
lMT_rV1c = dataArray{:, 11};
rMT_lV1p = dataArray{:, 12};
rMT_rV1p = dataArray{:, 13};
rMT_lV1c = dataArray{:, 14};
rMT_rV1c = dataArray{:, 15};
rMT_lMT = dataArray{:, 16};
lFEF_lV1p = dataArray{:, 17};
lFEF_rV1p = dataArray{:, 18};
lFEF_lV1c = dataArray{:, 19};
lFEF_rV1c = dataArray{:, 20};
lFEF_lMT = dataArray{:, 21};
lFEF_rMT = dataArray{:, 22};
rFEF_lV1p = dataArray{:, 23};
rFEF_rV1p = dataArray{:, 24};
rFEF_lV1c = dataArray{:, 25};
rFEF_rV1c = dataArray{:, 26};
rFEF_lMT = dataArray{:, 27};
rFEF_rMT = dataArray{:, 28};
rFEF_lFEF = dataArray{:, 29};
lSPL_lV1p = dataArray{:, 30};
lSPL_rV1p = dataArray{:, 31};
lSPL_lV1c = dataArray{:, 32};
lSPL_rV1c = dataArray{:, 33};
lSPL_lMT = dataArray{:, 34};
lSPL_rMT = dataArray{:, 35};
lSPL_lFEF = dataArray{:, 36};
lSPL_rFEF = dataArray{:, 37};
rSPL_lV1p = dataArray{:, 38};
rSPL_rV1p = dataArray{:, 39};
rSPL_lV1c = dataArray{:, 40};
rSPL_rV1c = dataArray{:, 41};
rSPL_lMT = dataArray{:, 42};
rSPL_rMT = dataArray{:, 43};
rSPL_lFEF = dataArray{:, 44};
rSPL_rFEF = dataArray{:, 45};
rSPL_lSPL = dataArray{:, 46};
lFFA_lV1p = dataArray{:, 47};
lFFA_rV1p = dataArray{:, 48};
lFFA_lV1c = dataArray{:, 49};
lFFA_rV1c = dataArray{:, 50};
lFFA_lMT = dataArray{:, 51};
lFFA_rMT = dataArray{:, 52};
lFFA_lFEF = dataArray{:, 53};
lFFA_rFEF = dataArray{:, 54};
lFFA_lSPL = dataArray{:, 55};
lFFA_rSPL = dataArray{:, 56};
rFFA_lV1p = dataArray{:, 57};
rFFA_rV1p = dataArray{:, 58};
rFFA_lV1c = dataArray{:, 59};
rFFA_rV1c = dataArray{:, 60};
rFFA_lMT = dataArray{:, 61};
rFFA_rMT = dataArray{:, 62};
rFFA_lFEF = dataArray{:, 63};
rFFA_rFEF = dataArray{:, 64};
rFFA_lSPL = dataArray{:, 65};
rFFA_rSPL = dataArray{:, 66};
rFFA_lFFA = dataArray{:, 67};
lOFA_lV1p = dataArray{:, 68};
lOFA_rV1p = dataArray{:, 69};
lOFA_lV1c = dataArray{:, 70};
lOFA_rV1c = dataArray{:, 71};
lOFA_lMT = dataArray{:, 72};
lOFA_rMT = dataArray{:, 73};
lOFA_lFEF = dataArray{:, 74};
lOFA_rFEF = dataArray{:, 75};
lOFA_lSPL = dataArray{:, 76};
lOFA_rSPL = dataArray{:, 77};
lOFA_lFFA = dataArray{:, 78};
lOFA_rFFA = dataArray{:, 79};
rOFA_lV1p = dataArray{:, 80};
rOFA_rV1p = dataArray{:, 81};
rOFA_lV1c = dataArray{:, 82};
rOFA_rV1c = dataArray{:, 83};
rOFA_lMT = dataArray{:, 84};
rOFA_rMT = dataArray{:, 85};
rOFA_lFEF = dataArray{:, 86};
rOFA_rFEF = dataArray{:, 87};
rOFA_lSPL = dataArray{:, 88};
rOFA_rSPL = dataArray{:, 89};
rOFA_lFFA = dataArray{:, 90};
rOFA_rFFA = dataArray{:, 91};
rOFA_lOFA = dataArray{:, 92};
lATC_lV1p = dataArray{:, 93};
lATC_rV1p = dataArray{:, 94};
lATC_lV1c = dataArray{:, 95};
lATC_rV1c = dataArray{:, 96};
lATC_lMT = dataArray{:, 97};
lATC_rMT = dataArray{:, 98};
lATC_lFEF = dataArray{:, 99};
lATC_rFEF = dataArray{:, 100};
lATC_lSPL = dataArray{:, 101};
lATC_rSPL = dataArray{:, 102};
lATC_lFFA = dataArray{:, 103};
lATC_rFFA = dataArray{:, 104};
lATC_lOFA = dataArray{:, 105};
lATC_rOFA = dataArray{:, 106};
rATC_lV1p = dataArray{:, 107};
rATC_rV1p = dataArray{:, 108};
rATC_lV1c = dataArray{:, 109};
rATC_rV1c = dataArray{:, 110};
rATC_lMT = dataArray{:, 111};
rATC_rMT = dataArray{:, 112};
rATC_lFEF = dataArray{:, 113};
rATC_rFEF = dataArray{:, 114};
rATC_lSPL = dataArray{:, 115};
rATC_rSPL = dataArray{:, 116};
rATC_lFFA = dataArray{:, 117};
rATC_rFFA = dataArray{:, 118};
rATC_lOFA = dataArray{:, 119};
rATC_rOFA = dataArray{:, 120};
rATC_lATC = dataArray{:, 121};
lV1_lV1p = dataArray{:, 122};
lV1_rV1p = dataArray{:, 123};
lV1_lV1c = dataArray{:, 124};
lV1_rV1c = dataArray{:, 125};
lV1_lMT = dataArray{:, 126};
lV1_rMT = dataArray{:, 127};
lV1_lFEF = dataArray{:, 128};
lV1_rFEF = dataArray{:, 129};
lV1_lSPL = dataArray{:, 130};
lV1_rSPL = dataArray{:, 131};
lV1_lFFA = dataArray{:, 132};
lV1_rFFA = dataArray{:, 133};
lV1_lOFA = dataArray{:, 134};
lV1_rOFA = dataArray{:, 135};
lV1_lATC = dataArray{:, 136};
lV1_rATC = dataArray{:, 137};
rV1_lV1p = dataArray{:, 138};
rV1_rV1p = dataArray{:, 139};
rV1_lV1c = dataArray{:, 140};
rV1_rV1c = dataArray{:, 141};
rV1_lMT = dataArray{:, 142};
rV1_rMT = dataArray{:, 143};
rV1_lFEF = dataArray{:, 144};
rV1_rFEF = dataArray{:, 145};
rV1_lSPL = dataArray{:, 146};
rV1_rSPL = dataArray{:, 147};
rV1_lFFA = dataArray{:, 148};
rV1_rFFA = dataArray{:, 149};
rV1_lOFA = dataArray{:, 150};
rV1_rOFA = dataArray{:, 151};
rV1_lATC = dataArray{:, 152};
rV1_rATC = dataArray{:, 153};
rV1_lV1 = dataArray{:, 154};
lmPFC_lV1p = dataArray{:, 155};
lmPFC_rV1p = dataArray{:, 156};
lmPFC_lV1c = dataArray{:, 157};
lmPFC_rV1c = dataArray{:, 158};
lmPFC_lMT = dataArray{:, 159};
lmPFC_rMT = dataArray{:, 160};
lmPFC_lFEF = dataArray{:, 161};
lmPFC_rFEF = dataArray{:, 162};
lmPFC_lSPL = dataArray{:, 163};
lmPFC_rSPL = dataArray{:, 164};
lmPFC_lFFA = dataArray{:, 165};
lmPFC_rFFA = dataArray{:, 166};
lmPFC_lOFA = dataArray{:, 167};
lmPFC_rOFA = dataArray{:, 168};
lmPFC_lATC = dataArray{:, 169};
lmPFC_rATC = dataArray{:, 170};
lmPFC_lV1 = dataArray{:, 171};
lmPFC_rV1 = dataArray{:, 172};
rmPFC_lV1p = dataArray{:, 173};
rmPFC_rV1p = dataArray{:, 174};
rmPFC_lV1c = dataArray{:, 175};
rmPFC_rV1c = dataArray{:, 176};
rmPFC_lMT = dataArray{:, 177};
rmPFC_rMT = dataArray{:, 178};
rmPFC_lFEF = dataArray{:, 179};
rmPFC_rFEF = dataArray{:, 180};
rmPFC_lSPL = dataArray{:, 181};
rmPFC_rSPL = dataArray{:, 182};
rmPFC_lFFA = dataArray{:, 183};
rmPFC_rFFA = dataArray{:, 184};
rmPFC_lOFA = dataArray{:, 185};
rmPFC_rOFA = dataArray{:, 186};
rmPFC_lATC = dataArray{:, 187};
rmPFC_rATC = dataArray{:, 188};
rmPFC_lV1 = dataArray{:, 189};
rmPFC_rV1 = dataArray{:, 190};
rmPFC_lmPFC = dataArray{:, 191};
lAG_lV1p = dataArray{:, 192};
lAG_rV1p = dataArray{:, 193};
lAG_lV1c = dataArray{:, 194};
lAG_rV1c = dataArray{:, 195};
lAG_lMT = dataArray{:, 196};
lAG_rMT = dataArray{:, 197};
lAG_lFEF = dataArray{:, 198};
lAG_rFEF = dataArray{:, 199};
lAG_lSPL = dataArray{:, 200};
lAG_rSPL = dataArray{:, 201};
lAG_lFFA = dataArray{:, 202};
lAG_rFFA = dataArray{:, 203};
lAG_lOFA = dataArray{:, 204};
lAG_rOFA = dataArray{:, 205};
lAG_lATC = dataArray{:, 206};
lAG_rATC = dataArray{:, 207};
lAG_lV1 = dataArray{:, 208};
lAG_rV1 = dataArray{:, 209};
lAG_lmPFC = dataArray{:, 210};
lAG_rmPFC = dataArray{:, 211};
rAG_lV1p = dataArray{:, 212};
rAG_rV1p = dataArray{:, 213};
rAG_lV1c = dataArray{:, 214};
rAG_rV1c = dataArray{:, 215};
rAG_lMT = dataArray{:, 216};
rAG_rMT = dataArray{:, 217};
rAG_lFEF = dataArray{:, 218};
rAG_rFEF = dataArray{:, 219};
rAG_lSPL = dataArray{:, 220};
rAG_rSPL = dataArray{:, 221};
rAG_lFFA = dataArray{:, 222};
rAG_rFFA = dataArray{:, 223};
rAG_lOFA = dataArray{:, 224};
rAG_rOFA = dataArray{:, 225};
rAG_lATC = dataArray{:, 226};
rAG_rATC = dataArray{:, 227};
rAG_lV1 = dataArray{:, 228};
rAG_rV1 = dataArray{:, 229};
rAG_lmPFC = dataArray{:, 230};
rAG_rmPFC = dataArray{:, 231};
rAG_lAG = dataArray{:, 232};
lMTG_lV1p = dataArray{:, 233};
lMTG_rV1p = dataArray{:, 234};
lMTG_lV1c = dataArray{:, 235};
lMTG_rV1c = dataArray{:, 236};
lMTG_lMT = dataArray{:, 237};
lMTG_rMT = dataArray{:, 238};
lMTG_lFEF = dataArray{:, 239};
lMTG_rFEF = dataArray{:, 240};
lMTG_lSPL = dataArray{:, 241};
lMTG_rSPL = dataArray{:, 242};
lMTG_lFFA = dataArray{:, 243};
lMTG_rFFA = dataArray{:, 244};
lMTG_lOFA = dataArray{:, 245};
lMTG_rOFA = dataArray{:, 246};
lMTG_lATC = dataArray{:, 247};
lMTG_rATC = dataArray{:, 248};
lMTG_lV1 = dataArray{:, 249};
lMTG_rV1 = dataArray{:, 250};
lMTG_lmPFC = dataArray{:, 251};
lMTG_rmPFC = dataArray{:, 252};
lMTG_lAG = dataArray{:, 253};
lMTG_rAG = dataArray{:, 254};
rMTG_lV1p = dataArray{:, 255};
rMTG_rV1p = dataArray{:, 256};
rMTG_lV1c = dataArray{:, 257};
rMTG_rV1c = dataArray{:, 258};
rMTG_lMT = dataArray{:, 259};
rMTG_rMT = dataArray{:, 260};
rMTG_lFEF = dataArray{:, 261};
rMTG_rFEF = dataArray{:, 262};
rMTG_lSPL = dataArray{:, 263};
rMTG_rSPL = dataArray{:, 264};
rMTG_lFFA = dataArray{:, 265};
rMTG_rFFA = dataArray{:, 266};
rMTG_lOFA = dataArray{:, 267};
rMTG_rOFA = dataArray{:, 268};
rMTG_lATC = dataArray{:, 269};
rMTG_rATC = dataArray{:, 270};
rMTG_lV1 = dataArray{:, 271};
rMTG_rV1 = dataArray{:, 272};
rMTG_lmPFC = dataArray{:, 273};
rMTG_rmPFC = dataArray{:, 274};
rMTG_lAG = dataArray{:, 275};
rMTG_rAG = dataArray{:, 276};
rMTG_lMTG = dataArray{:, 277};
lPCC_lV1p = dataArray{:, 278};
lPCC_rV1p = dataArray{:, 279};
lPCC_lV1c = dataArray{:, 280};
lPCC_rV1c = dataArray{:, 281};
lPCC_lMT = dataArray{:, 282};
lPCC_rMT = dataArray{:, 283};
lPCC_lFEF = dataArray{:, 284};
lPCC_rFEF = dataArray{:, 285};
lPCC_lSPL = dataArray{:, 286};
lPCC_rSPL = dataArray{:, 287};
lPCC_lFFA = dataArray{:, 288};
lPCC_rFFA = dataArray{:, 289};
lPCC_lOFA = dataArray{:, 290};
lPCC_rOFA = dataArray{:, 291};
lPCC_lATC = dataArray{:, 292};
lPCC_rATC = dataArray{:, 293};
lPCC_lV1 = dataArray{:, 294};
lPCC_rV1 = dataArray{:, 295};
lPCC_lmPFC = dataArray{:, 296};
lPCC_rmPFC = dataArray{:, 297};
lPCC_lAG = dataArray{:, 298};
lPCC_rAG = dataArray{:, 299};
lPCC_lMTG = dataArray{:, 300};
lPCC_rMTG = dataArray{:, 301};
rPCC_lV1p = dataArray{:, 302};
rPCC_rV1p = dataArray{:, 303};
rPCC_lV1c = dataArray{:, 304};
rPCC_rV1c = dataArray{:, 305};
rPCC_lMT = dataArray{:, 306};
rPCC_rMT = dataArray{:, 307};
rPCC_lFEF = dataArray{:, 308};
rPCC_rFEF = dataArray{:, 309};
rPCC_lSPL = dataArray{:, 310};
rPCC_rSPL = dataArray{:, 311};
rPCC_lFFA = dataArray{:, 312};
rPCC_rFFA = dataArray{:, 313};
rPCC_lOFA = dataArray{:, 314};
rPCC_rOFA = dataArray{:, 315};
rPCC_lATC = dataArray{:, 316};
rPCC_rATC = dataArray{:, 317};
rPCC_lV1 = dataArray{:, 318};
rPCC_rV1 = dataArray{:, 319};
rPCC_lmPFC = dataArray{:, 320};
rPCC_rmPFC = dataArray{:, 321};
rPCC_lAG = dataArray{:, 322};
rPCC_rAG = dataArray{:, 323};
rPCC_lMTG = dataArray{:, 324};
rPCC_rMTG = dataArray{:, 325};
rPCC_lPCC = dataArray{:, 326};

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;




%% Create correlation structures

vars=whos;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'_lMT$','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lMT_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lFEF','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lFEF_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lSPL','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lSPL_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'_rMT$','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    elseif (~isempty(regexp(vars(i).name,'rMT_lMT$','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rMT_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rFEF','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rFEF_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rSPL','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rSPL_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rFFA','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rFFA_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rOFA','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rOFA_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rV1','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rV1_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rATC','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rATC_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lFFA','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lFFA_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lOFA','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lOFA_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lV1','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lV1_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lATC','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lATC_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lmPFC','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lmPFC_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rmPFC','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rmPFC_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lAG','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lAG_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rAG','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rAG_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lMTG','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lMTG_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rMTG','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rMTG_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lPCC','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lPCC_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rPCC','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rPCC_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lV1c','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lV1c_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'lV1p','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
lV1p_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rV1c','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rV1c_corrs_YC = temp_struct;

temp_struct = struct;
for i=1:size(vars,1)
    if(~isempty(regexp(vars(i).name,'rV1p','match')))
        temp_struct.(vars(i).name)=eval(vars(i).name);
    end
end
rV1p_corrs_YC = temp_struct;

%% Save corrs to file

save('indiv_subj_corrs_z_YC', '*corrs_YC')
