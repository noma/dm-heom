// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/constants.hpp"

namespace heom {
namespace constants {

const long double eta_pade[] = {
	2.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.967175818975844784692668313888340351402, 
	1.032824181024155215307331686111659648598, 0, 0, 0, 0, 0, 0, 0, 0, 
	11.19885869716451815284926959121234309769, 1.300914314466752848423178399066599315543, 
	1.000226988368728998727552009721057586765, 0, 0, 0, 0, 0, 0, 0, 
	18.0790807744036381353279225690526377613, 1.905605223763320862837482220696034638651, 
	1.015313588078062122619093914564704203251, 1.000000413754978879215501295686623396803, 0, 0, 0, 0, 
	0, 0, 26.58655772186703141593168716080327785988, 2.800146619916205478116315260146686752654, 
	1.11303340660549395961527574088616034253, 1.000262251326821823459883045067529050956, 
	1.000000000284447322876838793096345993981, 0, 0, 0, 0, 0, 
	36.71706591275454016953866376131015396352, 3.915524608369778138857104011713377210806, 
	1.36060762304475364762294649920250800437, 1.006800099494035677225303524345587578974, 
	1.000001756336801785832037797152047562438, 1.000000000000090580923944406276325679894, 0, 0, 0, 0, 
	48.46945901876815939911017902590799133678, 5.221358700226750582830885303574125272024, 
	1.76101413086885578701660606207054290627, 1.048001360645944284478966960182497444605, 
	1.000166783806547433957971903891016917251, 1.000000005683742497388628469835034317229, 
	1.000000000000000015216762274538791805838, 0, 0, 0, 61.84333899016684164883279780411421181762, 
	6.709368642028293738757159220036413607145, 2.276634739389544433995764512332145784496, 
	1.167662545543253700136635145976136601895, 1.002993112803799185748365319635086694072, 
	1.000001970058157018499040773410845478966, 1.000000000010110274028760159309306116011, 
	1.000000000000000000001477065185853899797, 0, 0, 76.83854093567190329220602502323513693626, 
	8.377617889976198321037045566932019912232, 2.875743265352905559741498484957445572073, 
	1.386916450076667643674849716161379189731, 1.021092605957815260732482889860564891757, 
	1.000088840018504469980284336164406834212, 1.000000012945994679324047582718370173438, 
	1.000000000000010773303766311318241171588, 1.000000000000000000000000088652435318708, 0, 
	93.45498712043917746792411329602543487242, 10.22570574466955029693172438665880309108, 
	3.544527893692194584758656451673749774362, 1.692413746219203914341865930002582507947, 
	1.081051704862394372162728930732659549555, 1.00131233986343205708864037801373384889, 
	1.000001450202211357298203645556646479551, 1.000000000051835942172273790638144172787, 
	1.000000000000000007321793190694778206164, 1.000000000000000000000000000003467497243
};

const long double xi_pade[] = {
	7.745966692414833770358530799564799221666, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	19.49961875292266554276226912855528560143, 6.305939144224808474630860431551181929313, 0, 0, 0, 0, 
	0, 0, 0, 0, 36.11928942365857887916363930568047443027, 12.95828669566379113836433196980720441554, 
	6.283290335721270010976001332769855412497, 0, 0, 0, 0, 0, 0, 0, 
	57.78794000633463555484777318092706162218, 20.56259756752951982974170033395653675754, 
	12.57995038432843325319503152817888327099, 6.283185452295109512774977729104135871222, 0, 0, 0, 0, 
	0, 0, 84.53692559026397942121263237089439381387, 29.57927623718706979609581199085607351807, 
	19.00468956251068382661480700055507075954, 12.56654232705146493528755058935472519371, 
	6.28318530726039477349210954220860566044, 0, 0, 0, 0, 0, 
	116.3739347608330350073189696097089226231, 40.21293928856715127416984659035932828755, 
	25.82783232953028238453829864275777795836, 18.8562400267953494905521831656114592347, 
	12.56637153917480629098401692212665860923, 6.28318530717960815500107855468522104325, 0, 0, 0, 0, 
	153.3015628393653726331899655194706876159, 52.52412157944826796685914142534634337648, 
	33.34898703451500336730287935111300571215, 25.19900528903505310331751346668259596763, 
	18.84968412056844144448843441419311682201, 12.56637061688110415733890859371956448695, 
	6.283185307179586480075153567585578254575, 0, 0, 0, 195.3208905359343070977627172198676228279, 
	66.52807871550278259674646939491255055216, 41.76755706259573674247409323452395574559, 
	31.72895561746715247842724612158554657776, 25.13585329999435832854980649807171916496, 
	18.84955718260564459859203580283272390792, 12.56637061436306356083304954158209077981, 
	6.283185307179586476925556323382874292731, 0, 0, 242.43243633332578223865362314113358865, 
	82.22872388545609111729841286466997903147, 51.17127400682872814625269028495752604681, 
	38.64346654893615544641179566060655896253, 31.44483562147124840697966631064127684252, 
	25.13281634168110515436812140082770700199, 18.84955592869136612218434285960431278578, 
	12.56637061435917662145272429955102313797, 6.283185307179586476925286781031145512608, 0, 
	294.6364754696112018848669150774790893224, 99.62709625387022543419284510745749131322, 
	61.58845079917599582892471747985579618666, 46.1203736291541763005458261878858666512, 
	37.84390364517245088090747934823213310617, 31.41734004510526384373874233227272273452, 
	25.13274227487208277473085335076650873783, 18.84955592156405141516275123243075673053, 
	12.56637061435917295608825501073938367082, 6.283185307179586476925286766559517951084
};

} // namespace constants
} // namespace heom
