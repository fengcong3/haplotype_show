function getData() {
    let str = "name,<10,10-19,20-29,30-39,40-49,50-59,60-69,70-79,≥80\nAL,598478,638789,661666,603013,625599,673864,548376,316598,174781\nAK,106741,99926,120674,102008,91539,104569,70473,28422,12503\nAZ,892083,912735,939804,857054,833290,834858,737884,466153,254716\nAR,392177,397185,399698,372998,370157,395070,329734,197985,113468\nCA,5038433,5170341,5809455,5354112,5179258,5042094,3737461,2011678,1311374\nCO,690830,697447,780508,766382,705450,725661,563376,274466,155175\nCT,399369,481065,462323,424890,496265,546361,400995,217827,159475\nDE,112177,117854,127554,114063,117588,133331,110822,65369,35937\nDC,74377,62783,136976,121520,80570,74779,56984,31362,19658\nFL,2211012,2331102,2597830,2416176,2575576,2762983,2404659,1615547,1019566\nGA,1363631,1421557,1418696,1357210,1404698,1337985,998253,528108,269182\nHI,176484,163559,204336,187590,176904,188438,164957,85345,66060\nID,236658,239509,218684,209500,194678,205170,179429,97621,54234\nIL,1619682,1715984,1789739,1721954,1697069,1773366,1326121,728821,478948\nIN,863029,905594,905590,827086,844059,911778,704523,384788,243131\nIA,401712,418667,419456,383638,370719,427554,344037,197223,143583\nKS,401751,402092,406956,368732,344427,389834,300759,166104,117637\nKY,555615,575866,593819,558201,580553,623164,495736,273961,155074\nLA,622061,613633,683606,615411,571991,631936,488846,266123,152063\nME,137954,155774,156359,147695,176908,215787,179540,97899,62007\nMD,741952,764730,815346,784097,815875,862778,636309,330736,208079\nMA,737748,862371,971340,847306,916106,979128,737805,401931,288408\nMI,1181424,1324071,1338179,1162186,1283122,1454462,1148131,619722,398303\nMN,711604,714399,728222,715583,692201,782655,577313,312906,215985\nMS,400288,421329,414195,374724,377165,400164,319443,181195,100689\nMO,763948,792935,831725,763002,750989,857534,668878,388086,242554\nMT,126020,126294,136346,125004,116502,149800,130977,70528,41920\nNE,263518,257610,260646,244236,222479,250911,195705,107650,78504\nNV,369362,360263,392834,390261,387272,373757,309651,173499,82273\nNH,138762,167495,167554,151409,182703,217950,164287,84791,52552\nNJ,1079136,1153625,1139927,1143452,1254602,1307263,946399,523620,367432\nNM,276468,282662,289801,260579,244346,280363,239044,135013,74393\nNY,2319945,2445591,2894266,2605355,2617327,2755620,2095207,1160055,804091\nNC,1250298,1310398,1350242,1268976,1357746,1356117,1095320,609234,342497\nND,99477,91069,124509,94713,80327,98688,73825,41348,32206\nOH,1422838,1530264,1535538,1398724,1490959,1677794,1320776,728158,481890\nOK,534481,522282,552528,501392,469410,512850,404704,239887,138055\nOR,474456,485345,538596,537767,507826,534421,490894,255809,157153\nPA,1458931,1608018,1712448,1520409,1645291,1881378,1491536,850897,615069\nRI,111377,136885,153674,126503,137892,156127,117653,63359,51021\nSC,599591,619144,667523,596491,619792,663408,579856,322073,166727\nSD,120366,113383,116748,105499,96288,117012,92824,50398,38540\nTN,818404,842873,895337,837313,866343,904272,741045,414939,227483\nTX,3983091,3910528,3946447,3770534,3545746,3344930,2431494,1291486,732179\nUT,513515,479126,465219,436010,328569,301596,230007,123674,70711\nVT,63768,79351,81765,70092,79982,99521,82136,42978,26656\nVA,1033629,1065461,1170634,1112111,1134928,1162028,881763,475141,274606\nWA,895790,882812,1004428,970613,921205,970407,784208,401094,242589\nWV,207017,218547,232027,220494,238218,269346,243108,138134,79201\nWI,705648,755224,760961,714479,732280,848672,645015,350772,241747\nWY,78217,75535,82898,76912,68464,81018,67484,32819,19682\nPR,389428,479749,480184,441842,456009,452503,411924,268783,148963"
    return {
        'data': str
    }
}

function add(a, b) {
    return a + b;
}



var Utils = {
    getX() {
        return 12;
    },
    getY() {
        return 46;
    }
}

function constant$a(x) {
    return function constant() {
        return x;
    };
}

function none$1(series, order) {
    if (!((n = series.length) > 1)) return;
    for (var i = 1, j, s0, s1 = series[order[0]], n, m = s1.length; i < n; ++i) {
        s0 = s1, s1 = series[order[i]];
        for (j = 0; j < m; ++j) {
            s1[j][1] += s1[j][0] = isNaN(s0[j][1]) ? s0[j][0] : s0[j][1];
        }
    }
}

function none$2(series) {
    var n = series.length,
        o = new Array(n);
    while (--n >= 0) o[n] = n;
    return o;
}

function stackValue(d, key) {
    return d[key];
}

function stackSeries(key) {
    const series = [];
    series.key = key;
    return series;
}

function array$5(x) {
    return typeof x === "object" && "length" in x ?
        x // Array, TypedArray, NodeList, array-like
        :
        Array.from(x); // Map, Set, iterable, string, or anything else
}

function stack() {
    var keys = constant$a([]),
        order = none$2,
        offset = none$1,
        value = stackValue;

    function stack(data) {
        var sz = Array.from(keys.apply(this, arguments), stackSeries),
            i, n = sz.length,
            j = -1,
            oz;

        for (const d of data) {
            for (i = 0, ++j; i < n; ++i) {
                (sz[i][j] = [0, +value(d, sz[i].key, j, data)]).data = d;
            }
        }

        for (i = 0, oz = array$5(order(sz)); i < n; ++i) {
            sz[oz[i]].index = i;
        }

        offset(sz, oz);
        return sz;
    }

    stack.keys = function (_) {
        return arguments.length ? (keys = typeof _ === "function" ? _ : constant$a(Array.from(_)), stack) : keys;
    };

    stack.value = function (_) {
        return arguments.length ? (value = typeof _ === "function" ? _ : constant$a(+_), stack) : value;
    };

    stack.order = function (_) {
        return arguments.length ? (order = _ == null ? none$2 : typeof _ === "function" ? _ : constant$a(Array.from(_)), stack) : order;
    };

    stack.offset = function (_) {
        return arguments.length ? (offset = _ == null ? none$1 : _, stack) : offset;
    };

    return stack;
}

function expand(series, order) {
    if (!((n = series.length) > 0)) return;
    for (var i, n, j = 0, m = series[0].length, y; j < m; ++j) {
        for (y = i = 0; i < n; ++i) y += series[i][j][1] || 0;
        if (y)
            for (i = 0; i < n; ++i) series[i][j][1] /= y;
    }
    none$1(series, order);
}

function countKeys(data) {
    return Object.keys(data)
}