/* Copyright (c) <2009> <Newton Game Dynamics>
* 
* This software is provided 'as-is', without any express or implied
* warranty. In no event will the authors be held liable for any damages
* arising from the use of this software.
* 
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely
*/

#include "dContainersStdAfx.h"
#include "dCRC.h"


static dCRCTYPE randBits0[] =
{
	7266447313870364031uLL,  4946485549665804864uLL, 16945909448695747420uLL, 16394063075524226720uLL,  
	4873882236456199058uLL, 14877448043947020171uLL,  6740343660852211943uLL, 13857871200353263164uLL,
	5249110015610582907uLL, 10205081126064480383uLL,  1235879089597390050uLL, 17320312680810499042uLL,
	16489141110565194782uLL,  8942268601720066061uLL, 13520575722002588570uLL, 14226945236717732373uLL, 

	9383926873555417063uLL, 15690281668532552105uLL, 11510704754157191257uLL, 15864264574919463609uLL, 
	6489677788245343319uLL,  5112602299894754389uLL, 10828930062652518694uLL, 15942305434158995996uLL,
	15445717675088218264uLL, 4764500002345775851uLL, 14673753115101942098uLL,  236502320419669032uLL, 
	13670483975188204088uLL, 14931360615268175698uLL, 8904234204977263924uLL, 12836915408046564963uLL, 

	12120302420213647524uLL, 15755110976537356441uLL,  5405758943702519480uLL, 10951858968426898805uLL, 
	17251681303478610375uLL,  4144140664012008120uLL, 18286145806977825275uLL, 13075804672185204371uLL, 
	10831805955733617705uLL,  6172975950399619139uLL, 12837097014497293886uLL, 12903857913610213846uLL,
	560691676108914154uLL,    1074659097419704618uLL, 14266121283820281686uLL, 11696403736022963346uLL, 

	13383246710985227247uLL,  7132746073714321322uLL, 10608108217231874211uLL, 9027884570906061560uLL, 
	12893913769120703138uLL, 15675160838921962454uLL,  2511068401785704737uLL, 14483183001716371453uLL, 
	3774730664208216065uLL,  5083371700846102796uLL,  9583498264570933637uLL, 17119870085051257224uLL, 
	5217910858257235075uLL, 10612176809475689857uLL,  1924700483125896976uLL,  7171619684536160599uLL,


	10949279256701751503uLL, 15596196964072664893uLL, 14097948002655599357uLL, 615821766635933047uLL, 
	5636498760852923045uLL, 17618792803942051220uLL, 580805356741162327uLL,   425267967796817241uLL, 
	8381470634608387938uLL, 13212228678420887626uLL, 16993060308636741960uLL, 957923366004347591uLL, 
	6210242862396777185uLL,  1012818702180800310uLL, 15299383925974515757uLL, 17501832009465945633uLL,

	17453794942891241229uLL, 15807805462076484491uLL,  8407189590930420827uLL,   974125122787311712uLL,
	1861591264068118966uLL, 997568339582634050uLL, 18046771844467391493uLL, 17981867688435687790uLL, 
	3809841506498447207uLL,  9460108917638135678uLL, 16172980638639374310uLL,   958022432077424298uLL, 
	4393365126459778813uLL, 13408683141069553686uLL, 13900005529547645957uLL, 15773550354402817866uLL,

	16475327524349230602uLL,  6260298154874769264uLL, 12224576659776460914uLL,  6405294864092763507uLL, 
	7585484664713203306uLL,  5187641382818981381uLL, 12435998400285353380uLL, 13554353441017344755uLL,
	646091557254529188uLL, 11393747116974949255uLL, 16797249248413342857uLL, 15713519023537495495uLL,
	12823504709579858843uLL,  4738086532119935073uLL, 4429068783387643752uLL,  585582692562183870uLL,

	1048280754023674130uLL, 6788940719869959076uLL, 11670856244972073775uLL, 2488756775360218862uLL,
	2061695363573180185uLL,  6884655301895085032uLL, 3566345954323888697uLL, 12784319933059041817uLL,
	4772468691551857254uLL,  6864898938209826895uLL,  7198730565322227090uLL, 2452224231472687253uLL,
	13424792606032445807uLL, 10827695224855383989uLL, 11016608897122070904uLL, 14683280565151378358uLL,

	7077866519618824360uLL, 17487079941198422333uLL, 3956319990205097495uLL,  5804870313319323478uLL,
	8017203611194497730uLL,  3310931575584983808uLL,  5009341981771541845uLL, 11772020174577005930uLL,
	3537640779967351792uLL,  6801855569284252424uLL, 17687268231192623388uLL, 12968358613633237218uLL, 
	1429775571144180123uLL, 10427377732172208413uLL, 12155566091986788996uLL, 16465954421598296115uLL, 

	12710429690464359999uLL, 9547226351541565595uLL, 12156624891403410342uLL,  2985938688676214686uLL,
	18066917785985010959uLL,  5975570403614438776uLL, 11541343163022500560uLL, 11115388652389704592uLL,
	9499328389494710074uLL,  9247163036769651820uLL,  3688303938005101774uLL, 2210483654336887556uLL,
	15458161910089693228uLL,  6558785204455557683uLL,  1288373156735958118uLL, 18433986059948829624uLL,

	3435082195390932486uLL, 16822351800343061990uLL,  3120532877336962310uLL, 16681785111062885568uLL,  
	7835551710041302304uLL, 2612798015018627203uLL, 15083279177152657491uLL, 6591467229462292195uLL,
	10592706450534565444uLL,  7438147750787157163uLL, 323186165595851698uLL, 7444710627467609883uLL,
	8473714411329896576uLL,  2782675857700189492uLL,  3383567662400128329uLL, 3200233909833521327uLL,

	12897601280285604448uLL, 3612068790453735040uLL, 8324209243736219497uLL, 15789570356497723463uLL,
	1083312926512215996uLL,  4797349136059339390uLL, 5556729349871544986uLL, 18266943104929747076uLL,
	1620389818516182276uLL, 172225355691600141uLL,  3034352936522087096uLL,  1266779576738385285uLL,
	3906668377244742888uLL,  6961783143042492788uLL, 17159706887321247572uLL,  4676208075243319061uLL,

	10315634697142985816uLL, 13435140047933251189uLL, 716076639492622016uLL,	13847954035438697558uLL,
	7195811275139178570uLL, 10815312636510328870uLL,  6214164734784158515uLL, 16412194511839921544uLL, 
	3862249798930641332uLL,  1005482699535576005uLL,  4644542796609371301uLL, 17600091057367987283uLL, 
	4209958422564632034uLL, 5419285945389823940uLL, 11453701547564354601uLL,  9951588026679380114uLL,

	7425168333159839689uLL,  8436306210125134906uLL, 11216615872596820107uLL,  3681345096403933680uLL,
	5770016989916553752uLL, 11102855936150871733uLL, 11187980892339693935uLL, 396336430216428875uLL, 
	6384853777489155236uLL,  7551613839184151117uLL, 16527062023276943109uLL, 13429850429024956898uLL,
	9901753960477271766uLL,  9731501992702612259uLL,  5217575797614661659uLL, 10311708346636548706uLL,

	15111747519735330483uLL, 4353415295139137513uLL,  1845293119018433391uLL, 11952006873430493561uLL,
	3531972641585683893uLL, 16852246477648409827uLL,	15956854822143321380uLL, 12314609993579474774uLL,
	16763911684844598963uLL, 16392145690385382634uLL,  1545507136970403756uLL, 17771199061862790062uLL,
	12121348462972638971uLL, 12613068545148305776uLL,   954203144844315208uLL,  1257976447679270605uLL, 

	3664184785462160180uLL,  2747964788443845091uLL, 15895917007470512307uLL, 15552935765724302120uLL,
	16366915862261682626uLL, 8385468783684865323uLL, 10745343827145102946uLL, 2485742734157099909uLL, 
	916246281077683950uLL, 15214206653637466707uLL, 12895483149474345798uLL,  1079510114301747843uLL,
	10718876134480663664uLL,  1259990987526807294uLL,  8326303777037206221uLL, 14104661172014248293uLL,
};


/*
static dCRCTYPE randBits1[] =
{
	7266447313870364031uLL, 4946485549665804864uLL, 1694590944869574742uLL, 1639406307552422672uLL, 
	4873882236456199058uLL, 1487744804394702017uLL, 6740343660852211943uLL, 1385787120035326316uLL,
	5249110015610582907uLL, 1020508112606448038uLL, 1235879089597390050uLL, 1732031268081049904uLL,
	1648914111056519478uLL, 8942268601720066061uLL, 1352057572200258857uLL, 1422694523671773237uLL, 

	9383926873555417063uLL, 1569028166853255210uLL, 1151070475415719125uLL, 1586426457491946360uLL, 
	6489677788245343319uLL, 5112602299894754389uLL, 1082893006265251869uLL, 1594230543415899599uLL,
	1544571767508821826uLL, 4764500002345775851uLL, 1467375311510194209uLL, 2365023204196690325uLL, 
	1367048397518820408uLL, 1493136061526817569uLL, 8904234204977263924uLL, 1283691540804656496uLL, 

	1212030242021364752uLL, 1575511097653735644uLL, 5405758943702519480uLL, 1095185896842689880uLL, 
	1725168130347861037uLL, 4144140664012008120uLL, 1828614580697782527uLL, 1307580467218520437uLL, 
	1083180595573361770uLL, 6172975950399619139uLL, 1283709701449729388uLL, 1290385791361021384uLL,
	5606916761089141545uLL, 1074659097419704618uLL, 1426612128382028168uLL, 1169640373602296334uLL, 

	1338324671098522724uLL, 7132746073714321322uLL, 1060810821723187421uLL, 9027884570906061560uLL, 
	1289391376912070313uLL, 1567516083892196245uLL, 2511068401785704737uLL, 1448318300171637145uLL, 
	3774730664208216065uLL, 5083371700846102796uLL, 9583498264570933637uLL, 1711987008505125722uLL, 
	5217910858257235075uLL, 1061217680947568985uLL, 1924700483125896976uLL, 7171619684536160599uLL,

	1094927925670175150uLL, 1559619696407266489uLL, 1409794800265559935uLL, 6158217666359330479uLL, 
	5636498760852923045uLL, 1761879280394205122uLL, 5808053567411623279uLL, 4252679677968172418uLL, 
	8381470634608387938uLL, 1321222867842088762uLL, 1699306030863674196uLL, 9579233660043475913uLL, 
	6210242862396777185uLL, 1012818702180800310uLL, 1529938392597451575uLL, 1750183200946594563uLL,

	1745379494289124122uLL, 1580780546207648449uLL, 8407189590930420827uLL, 9741251227873117122uLL,
	1861591264068118966uLL, 9975683395826340507uLL, 1804677184446739149uLL, 1798186768843568779uLL, 
	3809841506498447207uLL, 9460108917638135678uLL, 1617298063863937431uLL, 9580224320774242987uLL, 
	4393365126459778813uLL, 1340868314106955368uLL, 1390000552954764595uLL, 1577355035440281786uLL,

	1647532752434923060uLL, 6260298154874769264uLL, 1222457665977646091uLL, 6405294864092763507uLL, 
	7585484664713203306uLL, 5187641382818981381uLL, 1243599840028535338uLL, 1355435344101734475uLL,
	6460915572545291889uLL, 1139374711697494925uLL, 1679724924841334285uLL, 1571351902353749549uLL,
	1282350470957985884uLL, 4738086532119935073uLL, 4429068783387643752uLL, 5855826925621838703uLL,

	1048280754023674130uLL, 6788940719869959076uLL, 1167085624497207377uLL, 2488756775360218862uLL,
	2061695363573180185uLL, 6884655301895085032uLL, 3566345954323888697uLL, 1278431993305904181uLL,
	4772468691551857254uLL, 6864898938209826895uLL, 7198730565322227090uLL, 2452224231472687253uLL,
	1342479260603244580uLL, 1082769522485538398uLL, 1101660889712207090uLL, 1468328056515137835uLL,

	7077866519618824360uLL, 1748707994119842233uLL, 3956319990205097495uLL, 5804870313319323478uLL,
	8017203611194497730uLL, 3310931575584983808uLL, 5009341981771541845uLL, 1177202017457700593uLL,
	3537640779967351792uLL, 6801855569284252424uLL, 1768726823119262338uLL, 1296835861363323721uLL, 
	1429775571144180123uLL, 1042737773217220841uLL, 1215556609198678899uLL, 1646595442159829611uLL, 

	1271042969046435999uLL, 9547226351541565595uLL, 1215662489140341034uLL, 2985938688676214686uLL,
	1806691778598501095uLL, 5975570403614438776uLL, 1154134316302250056uLL, 1111538865238970459uLL,
	9499328389494710074uLL, 9247163036769651820uLL, 3688303938005101774uLL, 2210483654336887556uLL,
	1545816191008969322uLL, 6558785204455557683uLL, 1288373156735958118uLL, 1843398605994882962uLL,

	3435082195390932486uLL, 1682235180034306199uLL, 3120532877336962310uLL, 1668178511106288556uLL, 
	7835551710041302304uLL, 2612798015018627203uLL, 1508327917715265749uLL, 6591467229462292195uLL,
	1059270645053456544uLL, 7438147750787157163uLL, 3231861655958516989uLL, 7444710627467609883uLL,
	8473714411329896576uLL, 2782675857700189492uLL, 3383567662400128329uLL, 3200233909833521327uLL,

	1289760128028560444uLL, 3612068790453735040uLL, 8324209243736219497uLL, 1578957035649772346uLL,
	1083312926512215996uLL, 4797349136059339390uLL, 5556729349871544986uLL, 1826694310492974707uLL,
	1620389818516182276uLL, 1722253556916001411uLL, 3034352936522087096uLL, 1266779576738385285uLL,
	3906668377244742888uLL, 6961783143042492788uLL, 1715970688732124757uLL, 4676208075243319061uLL,

	1031563469714298581uLL, 1343514004793325118uLL, 7160766394926220168uLL, 1384795403543869755uLL,
	7195811275139178570uLL, 1081531263651032887uLL, 6214164734784158515uLL, 1641219451183992154uLL, 
	3862249798930641332uLL, 1005482699535576005uLL, 4644542796609371301uLL, 1760009105736798728uLL, 
	4209958422564632034uLL, 5419285945389823940uLL, 1145370154756435460uLL, 9951588026679380114uLL,

	7425168333159839689uLL, 8436306210125134906uLL, 1121661587259682010uLL, 3681345096403933680uLL,
	5770016989916553752uLL, 1110285593615087173uLL, 1118798089233969393uLL, 3963364302164288756uLL, 
	6384853777489155236uLL, 7551613839184151117uLL, 1652706202327694310uLL, 1342985042902495689uLL,
	9901753960477271766uLL, 9731501992702612259uLL, 5217575797614661659uLL, 1031170834663654870uLL,

	1511174751973533048uLL, 4353415295139137513uLL, 1845293119018433393uLL, 1195200687343049356uLL,
	3531972641585683893uLL, 1685224647764840982uLL, 1595685482214332138uLL, 1231460999357947477uLL,
	1676391168484459896uLL, 1639214569038538263uLL, 1545507136970403756uLL, 1777119906186279006uLL,
	1212134846297263897uLL, 1261306854514830577uLL, 9542031448443152088uLL, 1257976447679270605uLL, 

	3664184785462160180uLL, 2747964788443845091uLL, 1589591700747051230uLL, 1555293576572430212uLL,
	1636691586226168262uLL, 8385468783684865323uLL, 1074534382714510294uLL, 2485742734157099909uLL, 
	9162462810776839507uLL, 1521420665363746670uLL, 1289548314947434579uLL, 1079510114301747843uLL,
	1071887613448066366uLL, 1259990987526807294uLL, 8326303777037206221uLL, 1410466117201424829uLL,
};
*/


dCRCTYPE dCombineCRC (dCRCTYPE a, dCRCTYPE b)
{
	return (a << 8) ^ b;
}

// calculate a 32 bit crc of a string
dCRCTYPE dCRC64 (const char* const name, dCRCTYPE crcAcc)
{
	if (name) {
		const int bitshift = (sizeof (dCRCTYPE)<<3) - 8;
		for (int i = 0; name[i]; i ++) {
			char c = name[i];
			dCRCTYPE val = randBits0[((crcAcc >> bitshift) ^ c) & 0xff];
			crcAcc = (crcAcc << 8) ^ val;
		}
	}
	return crcAcc;
}


dCRCTYPE dCRC64 (const void* const buffer, int size, dCRCTYPE crcAcc)
{
	const unsigned char* const ptr = (unsigned char*)buffer;

	const int bitshift = (sizeof (dCRCTYPE)<<3) - 8;
	for (int i = 0; i < size; i ++) {
		char c = ptr[i];
		dCRCTYPE  val = randBits0[((crcAcc >> bitshift) ^ c) & 0xff];
		crcAcc = (crcAcc << 8) ^ val;
	}
	return crcAcc;
}



