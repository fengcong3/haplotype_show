// Special character escape
function htmlEscape(text) {
    return text.replace(/[<>"&]/g, function (match, pos, originalText) {
        switch (match) {
            case "<":
                return "&lt;";
            case ">":
                return "&gt;";
            case "&":
                return "&amp;";
            case "\"":
                return "&quot;";
        }
    });
}

// Date format
Date.prototype.format = function (fmt) {
    var o = {
        "M+": this.getMonth() + 1,                 //月份 
        "d+": this.getDate(),                    //日 
        "h+": this.getHours(),                   //小时 
        "m+": this.getMinutes(),                 //分 
        "s+": this.getSeconds(),                 //秒 
        "q+": Math.floor((this.getMonth() + 3) / 3), //季度 
        "S": this.getMilliseconds()             //毫秒 
    };
    if (/(y+)/.test(fmt)) {
        fmt = fmt.replace(RegExp.$1, (this.getFullYear() + "").substr(4 - RegExp.$1.length));
    }
    for (var k in o) {
        if (new RegExp("(" + k + ")").test(fmt)) {
            fmt = fmt.replace(RegExp.$1, (RegExp.$1.length == 1) ? (o[k]) : (("00" + o[k]).substr(("" + o[k]).length)));
        }
    }
    return fmt;
}




// Processing original data 
function handleResData(data) {
    let temp = [];
    let geneStructure = data["gene_structure"];
    let sampleNames = data["sample_name"];
    let variation = data["variation"];
    let exon = geneStructure["exon"];
    let upstream = geneStructure["upstream"]["0"];
    let ori = geneStructure["ori"];
    let gene = geneStructure["gene"]["0"];
    let geneName = geneStructure["gene_name"];
    let p5, p3;
    //判断5p_UTR和3p_UTR是否存在
    if (geneStructure["3p_UTR"].length === 0) {
        p3 = [];
    } else {
        p3 = geneStructure["3p_UTR"]["0"];
    }
    if (geneStructure["5p_UTR"].length === 0) {
        p5 = [];
    } else {
        p5 = geneStructure["5p_UTR"]["0"];
    }

    let gs_obj = {};
    if (ori === '+') {
        gs_obj["start"] = upstream["0"];
        gs_obj["end"] = gene["1"];
        gs_obj["geneEnd"] = gene["1"];
    } else if (ori === '-') {
        gs_obj["start"] = gene["0"];
        gs_obj["end"] = upstream["1"];
        gs_obj["geneEnd"] = gene["1"];
    }

    if (+(gs_obj["start"]) > 0) {
        if (p3.length != 0) {
            gs_obj["3p_UTR"] = p3["1"] - p3["0"] + 1 + gs_obj["start"];
            gs_obj["3p_UTR_start"] = p3["0"];
        } else {
            gs_obj["3p_UTR"] = 0;
            gs_obj["3p_UTR_start"] = 0;
        }
        if (p5.length != 0) {
            gs_obj["5p_UTR"] = p5["1"] - p5["0"] + 1 + gs_obj["start"];
            gs_obj["5p_UTR_start"] = p5["0"];
        } else {
            gs_obj["5p_UTR"] = 0;
            gs_obj["5p_UTR_start"] = 0;
        }
        gs_obj["upstream"] = upstream["1"] - upstream["0"] + 1 + gs_obj["start"];
        let _exon = [];
        exon.forEach(d => {
            _exon.push([d["0"], d["1"] + gs_obj["start"]]);

        })
        gs_obj["exon"] = _exon;
    } else {
        if (p3.length != 0) {
            gs_obj["3p_UTR"] = p3["1"] - p3["0"] + 1;
            gs_obj["3p_UTR_start"] = p3["0"];
        } else {
            gs_obj["3p_UTR"] = 0;
            gs_obj["3p_UTR_start"] = 0;
        }
        if (p5.length != 0) {
            gs_obj["5p_UTR"] = p5["1"] - p5["0"] + 1;
            gs_obj["5p_UTR_start"] = p5["0"];
        } else {
            gs_obj["5p_UTR"] = 0;
            gs_obj["5p_UTR_start"] = 0;
        }
        gs_obj["upstream"] = upstream["1"] - upstream["0"] + 1;
        gs_obj["exon"] = exon;
    }
    gs_obj["name"] = sampleNames;
    gs_obj["variation"] = variation;
    gs_obj["ori"] = ori;
    gs_obj["geneName"] = geneName;
    temp.push(gs_obj);
    return temp;
}

// download data 
function handleClickDownloadData() {
    let cls = ".download";
    let getGeneName = $(".showGeneName").text()
    let data = {
        "name": getGeneName,
        "act": "down"
    };
    console.log(getGeneName);
    // http://hap.bioinf.club/
    if (getGeneName != "") {
        $.ajax({
            type: "POST",
            url: "http://hap.bioinf.club/cgi-bin/down.py",
            data: data,
            contentType: "application/x-www-form-urlencoded",
            header: {
                "Access-Control-Allow-Origin": "*"
            },
            async: true,
            beforeSend: function () {
                setButtonState(cls);
            },
            success: function (data) {
                console.log(data);
                downloadFile(data, getGeneName + ".tar.gz");
                removeButtonState(cls);
            },
            error: function (err) {
                console.log(err);
                removeButtonState(cls);
            }
        });
    } else {
        return;
    }
}

function downloadFile(url, filename) {
    // 创建隐藏的可下载链接
    var eleLink = document.createElement('a');
    //设置文件下载名
    eleLink.download = filename;
    //隐藏元素
    eleLink.style.display = 'none';
    // // 字符内容转变成blob地址
    // var blob = new Blob([url]);
    // eleLink.href = URL.createObjectURL(blob);
    eleLink.href = url;
    // // 触发点击
    document.body.appendChild(eleLink);
    eleLink.click();
    // // 然后移除
    document.body.removeChild(eleLink);
    // // 释放blob对象
    // URL.revokeObjectURL(url);
};

// click btn to open new window and show gene tree
function handleClickToDrawTree(cls, data) {
    let btnState = $(cls).attr("disabled");
    if (btnState === "disabled") {
        return;
    } else {
        let showGeneName = $(".showGeneName").text();
        console.log(showGeneName);
        $.ajax({
            type: "POST",
            url: "http://hap.bioinf.club/cgi-bin/toiTol.py",
            data: data,
            contentType: "application/x-www-form-urlencoded",
            async: true,
            beforeSend: function () {
                setButtonState(cls);
            },
            success: function (data) {
                console.log(data);
                window.open(data);
                removeButtonState(cls);
            },
            error: function (err) {
                console.log(err);
                removeButtonState(cls);
            }
        });
    }
}

function handleClickToDrawGeneTree() {
    let cls = ".gene-phy";
    let showGeneName = $(".showGeneName").text()
    let data = {
        "name": showGeneName,
        "act": "tree1"
    };
    handleClickToDrawTree(cls, data);
}

function handleClickToDrawHapTree() {
    let cls = ".hap-phy";
    let showGeneName = $(".showGeneName").text();
    let data = {
        "name": showGeneName,
        "act": "tree2"
    };
    handleClickToDrawTree(cls, data);
}

// show loading effect
function loadingEffect() {
    var loading = $('#loading-effect');
    //loading.hide();
    $(document).ajaxStart(function () {
        loading.show();
    }).ajaxStop(function () {
        loading.hide();
    });
}

// Set button state when making Ajax request
function setButtonState(cls) {
    $(cls).attr("disabled", true);
    $(cls).css("opacity", 0.6);
}

// remove button state when Ajax request done
function removeButtonState(cls) {
    $(cls).removeAttr("disabled");
    $(cls).css("opacity", 1);
}

// show table 
function csvDataToHtml(tab) {
    let geneName = $(".showGeneName").text()
    let csvPath = "../gene_data/" + geneName + "." + tab + ".csv";
    // let csvPath = "../gene_data/TraesCS2D02G079600." + tab + ".csv"
    if (tab === "snp") {
        d3.csv(csvPath).then(function (data) {
            showTable(data);
        })
    } else if (tab === "sv") {
        d3.csv(csvPath).then(function (data) {
            showTable(data);
        })
    } else if (tab === "indel") {
        d3.csv(csvPath).then(function (data) {
            showTable(data);
        })
    }
}

function showTable(data) {
    columns = data.columns;
    let theadHtml = "";
    let tbodyHtml = ""
    theadHtml += "<tr>";
    columns.forEach(d => {
        theadHtml += "<td>" + d + "</td>";
    })
    theadHtml += "</tr>";
    let len = columns.length;
    data.forEach(d => {
        count = 1;
        tbodyHtml += "<tr>";
        columns.forEach(di => {
            if (count === len) {
                tbodyHtml += "<td class='last-td'>" + htmlEscape(d[di]) + "</td>";
            } else {
                tbodyHtml += "<td>" + htmlEscape(d[di]) + "</td>";
            }
            count += 1;
        })
        tbodyHtml += "</tr>";
    })
    $(".tbl-head-data").html(theadHtml);
    $(".tbl-body-data").html(tbodyHtml);
    handleClickTableCell();
}

// show th td title
function showTableText() {
    $(".tbl-content").mouseover(function () {
        $("tbody tr").each(function (i) {
            $(this).children('td').each(function (j) {
                // 设置title属性
                $(this).attr('title', $(this).text());
            });
        })
    })
}


function handleClickTableCell () {
    $("td.last-td").each(function (i) {
        $(this).on("click", function () {
            let text = $(this).text();
            text = text.trim().replace(/\//g, " ").trim();
            $(".filter-samplename").val(text);
        })
    })
}

function saveDataToFile () {
    let text = $(".filter-samplename").val();
    let time = new Date().format("MM-dd-hh:mm:ss");
    let strFile = time + "-saveFilterFile.txt";
    let element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    element.setAttribute('download', strFile);
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
    $(".filter-samplename").val("");
}


// enter to get other gene
function getResData(name) {
    let _url = "../gene_data/" + name + ".json";
    let _data = null;
    $.ajax({
        type: "GET",
        url: _url,
        dataType: "json",
        contentType: "application/json; charset=utf-8",
        async: false,
        success: function (resData) {
            _data = resData;
        },
        error: function (err) {
            console.log(err);
        }
    });
    return _data;
}

function handleSubmitRes(data, name) {
    response = getResData(name);
    if (response != null) {
        dealData = handleResData(response);
        return dealData;
    } else {
        alert("data not exist!");
        return response;
    }
}

function changeShowInfo(flag) {
    if (flag) {
        $(".snp").each(function () {
            $(this).attr("opacity", 1)
        });
    } else {
        $(".snp").each(function () {
            $(this).attr("opacity", 0)
        });
    }
}

function showCurrentPageSNPMutationInfo() {
    $(function () {
        $(".switch").change(function () {
            let change = $(this).prop('checked');
            changeShowInfo(change);
        })
    });
}
