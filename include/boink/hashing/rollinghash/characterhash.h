#ifndef CHARACTERHASH
#define CHARACTERHASH

#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>


template <typename hashvaluetype> hashvaluetype maskfnc(int bits) {
  assert(bits > 0);
  assert(bits <= sizeof(hashvaluetype) * 8);
  hashvaluetype x = static_cast<hashvaluetype>(1) << (bits - 1);
  return x ^ (x - 1);
}

template <typename hashvaluetype = uint64_t, typename chartype = unsigned char>
class CharacterHash {
public:
  CharacterHash(hashvaluetype maxval) {}

  enum { nbrofchars = 1 << (sizeof(chartype) * 8) };

  // these are the random integers that need regenerating if you change
  // the number of bits (wordsize) in cyclichash.h
  const uint64_t hashvalues[256] = {
      2821057300106699797ull,  11325163037179593313ull, 5416224960711267644ull,
      7528795173630926414ull,  13801220329386839123ull, 9571312227480569745ull,
      13287641504033217005ull, 17026797478042468886ull, 12854026315321965989ull,
      6374470260133235413ull,  2204210189948696549ull,  16294121409648186942ull,
      52948147103884203ull,    6415147802595498673ull,  16198407629368318724ull,
      11766006283457710847ull, 5173278555550225963ull,  291850381695702383ull,
      9067924282960618557ull,  7697914929610629098ull,  12522305300226911288ull,
      9498027593538051549ull,  17257649845518499103ull, 14516692441855410225ull,
      8715626173811934569ull,  5325903972055096315ull,  16166526287333920336ull,
      2569195063296827853ull,  17109513377892528534ull, 13886494041707436679ull,
      14893740505833355047ull, 4234952225844646826ull,  14358747021030577240ull,
      1887737051715844916ull,  2330098280227001344ull,  9188490525853320411ull,
      6025534915685878506ull,  17424697768797008296ull, 7637984810305105866ull,
      12367942525066656103ull, 9840387827592398565ull,  10533248810674727780ull,
      13207324779243726126ull, 7315623160354210225ull,  11948867329934421237ull,
      3644170777857093101ull,  9676205714709334695ull,  14771130086372921028ull,
      2292394639782612013ull,  10197765460546535208ull, 14073391231741944415ull,
      8621873250296467544ull,  6590463839587788851ull,  13850038707066944250ull,
      18429044708838245451ull, 215258643983951745ull,   1884062495295395303ull,
      12663706866596696824ull, 14449601608250063354ull, 3435894449503459935ull,
      10870768431824442842ull, 11504721573538595258ull, 12186578141009955262ull,
      16198870037728999207ull, 10669405064869886793ull, 16664410744025174816ull,
      18397359403060890687ull, 15956807086001210932ull, 14617266033040342736ull,
      13852442956164108892ull, 10306008263506153018ull, 9404339731978646439ull,
      8655922979817266314ull,  12809422002881587418ull, 4261482358731270443ull,
      7732774012761317500ull,  16502580177155673326ull, 15611680917896763796ull,
      10780076607752328646ull, 1542572181126078483ull,  3654327363130598725ull,
      16794808258298155332ull, 12369059170961750015ull, 3132818342866713807ull,
      836480985777824379ull,   13906561201785132592ull, 9858196578848854067ull,
      4372375965812152942ull,  14254521215283273274ull, 7596145618544622863ull,
      16342747097867571709ull, 1980540366340112804ull,  8262177221342635774ull,
      13125315094409684063ull, 17203349866499262586ull, 16759515313626684424ull,
      15342020348994741163ull, 2706664420015103890ull,  11384314333735679458ull,
      15295060863465421772ull, 6664089646357471419ull,  14233415411242489719ull,
      12908263413292462998ull, 1353330147242825292ull,  15303480952098533508ull,
      337359223552501433ull,   1568788723399388210ull,  17644968499059613652ull,
      3847406208385860132ull,  17523265639259566024ull, 12238335531163888847ull,
      3903849605456502145ull,  8263775787236913822ull,  12389831137335901677ull,
      1292775400312032391ull,  3771478312966411459ull,  4541792594786962408ull,
      14808262216606926175ull, 2581919202420840491ull,  166168189331466451ull,
      6554319341342777275ull,  17145090797639572225ull, 12732027427133868456ull,
      14559633844336852506ull, 909502806480392211ull,   5778876322733837940ull,
      16859623521438723622ull, 5134076461472812216ull,  12354633801538016303ull,
      17850322572567648614ull, 5820236151567115665ull,  6576849339431981618ull,
      10590600560778871269ull, 10940127808752622216ull, 8565035553835162937ull,
      1774065049477478434ull,  4814212482303396695ull,  4275189062710062116ull,
      630349411251672937ull,   6436038361922976831ull,  16664238356166903433ull,
      2398603741284037446ull,  11468261837055649518ull, 17514985377544487208ull,
      3177974761867308872ull,  7692698080149068968ull,  13837720152036376942ull,
      5506266379279335520ull,  5308518106289494774ull,  16526183736085529801ull,
      357257327638821318ull,   1181833586544916844ull,  16664823884918697019ull,
      2529708085974818140ull,  18394812950618671631ull, 9885543901797367607ull,
      10673912057166833585ull, 522162561327397618ull,   8433940106549935649ull,
      2621653722353075094ull,  9722505992377914047ull,  1711872796729989241ull,
      17798157461229462354ull, 16913106389905575378ull, 7896888628990885297ull,
      15864715499398966611ull, 14079164878488738891ull, 3050246499307582612ull,
      1703346441458602658ull,  7310689803047170057ull,  11498964545764371542ull,
      17202646980119590472ull, 9680839911901944918ull,  7640258924910894043ull,
      6529023846032485156ull,  483507376575727765ull,   13079525949120166236ull,
      10170150048283242299ull, 2109960879517206ull,     16760628739394415400ull,
      10582518701122574096ull, 4355072925142789383ull,  5781655723244612819ull,
      7768064507395797426ull,  921478671905296912ull,   3179121294549908816ull,
      9975781372392684305ull,  549738065842820422ull,   12640068765622852389ull,
      8229728960599470880ull,  9375512254211572010ull,  12744798652893158085ull,
      2363478901711075627ull,  2116689541286966759ull,  8180262042010236006ull,
      4249443692261027106ull,  2707168392994699624ull,  5514287253458790481ull,
      1814178923525761966ull,  5577051504293005778ull,  12771098131382550559ull,
      12032790465375591411ull, 2535960658173991567ull,  10820030687998681802ull,
      4979290776778324590ull,  16457002048535804954ull, 12244850006062512789ull,
      7609998496390090313ull,  14104243468864345707ull, 720433609917385434ull,
      13392299073419520539ull, 4097474698105279986ull,  8299413979248706606ull,
      16290789562718246821ull, 7335842782667280191ull,  16257887162134443769ull,
      17861273412307426758ull, 13701592676399062223ull, 17373825138829494763ull,
      7944986297238901628ull,  9835163643181430505ull,  366723710453002021ull,
      15396126198530254709ull, 11511321399930066684ull, 3540964612902366865ull,
      5036883055268744511ull,  4962193520082228875ull,  14120707866588336117ull,
      18241267409277997530ull, 17669943049921343994ull, 9932182549666343663ull,
      9558227449175976386ull,  1092844097255405105ull,  1904183829110741896ull,
      7155690151231400523ull,  10439234361950944381ull, 1361222343450861792ull,
      9095344136183365091ull,  4898471236227651240ull,  2589681900905316607ull,
      12762880031044612488ull, 9939678776623939418ull,  505211974413939563ull,
      9972319229645889248ull,  7319067952382618431ull,  8681790484606533729ull,
      6965201783339143147ull,  7006134092273617374ull,  5342723917871999735ull,
      984365199659547641ull,   5559120030223996309ull,  9939414356931300058ull,
      1217488313221569951ull,  17323850086388522307ull, 15532728393800726547ull,
      13935833578031356880ull};
};

#endif
