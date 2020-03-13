function get_js_version() {
    return "JS2020-03-13";
}

function poll_progress() {
    var xmlhttp = new XMLHttpRequest();
    var url = 'progress/' + encodeURIComponent(uuid);

    xmlhttp.onreadystatechange = function () {
        if (xmlhttp.readyState == 4 && xmlhttp.status == 200) {
            var data = xmlhttp.response;

            try {
                process_progress_event(data);

                if (data.total == 0 || data.completed != data.total)
                    setTimeout(function () {
                        poll_progress();
                    }, 250);
            } catch (e) {
                console.log(e);
            }
        }

        if (xmlhttp.readyState == 4 && xmlhttp.status == 204) {
            console.log("data not found");
            $("#progress").remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    }

    xmlhttp.open("POST", url, true);
    xmlhttp.responseType = 'json';
    xmlhttp.timeout = 0;
    xmlhttp.send();
}

function poll_status() {
    var xmlhttp = new XMLHttpRequest();
    var url = 'status/' + encodeURIComponent(uuid);

    xmlhttp.onreadystatechange = function () {
        if (xmlhttp.readyState == 4 && xmlhttp.status == 200) {
            fetch_plots();// fetch_plots() or fetch_json_data()
        }

        if (xmlhttp.readyState == 4 && xmlhttp.status == 204) {
            console.log("data not found");
            $("#progress").remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }

        // repeat the poll until success
        if (xmlhttp.readyState == 4 && xmlhttp.status == 202) {
            setTimeout(function () {
                poll_status();
            }, 250);
        }
    }

    xmlhttp.open("POST", url, true);
    xmlhttp.timeout = 0;
    xmlhttp.send();
}

function process_progress_event(data) {
    if (data.type == "progress") {
        session_data = document.getElementById('session-data');
        va_count = parseInt(session_data.getAttribute('data-search-count'));
        //console.log("json:", data);
        $('#completed').html("completed: " + data.completed + "/" + va_count);

        var PROGRESS_VARIABLE = data.completed;
        var progress = PROGRESS_VARIABLE / va_count;
        var PROGRESS_INFO = "&nbsp;" + numeral(progress).format('0.0%');

        var speed = data.completed / data.elapsed;
        var remaining_time = (va_count - data.completed) / speed;
        //console.log("speed:", speed, "remaining:", remaining_time);
        if (remaining_time > 1)
            PROGRESS_INFO += ", est. time to completion: " + numeral(remaining_time).format('00:00:00');

        $("#progress-bar")
            .attr("aria-valuenow", PROGRESS_VARIABLE)
            .css("width", (100.0 * progress) + "%")
            .html(PROGRESS_INFO);

        if (PROGRESS_VARIABLE == va_count) {
            //close the websocket connection
            //gaia_ws.close();//keep the connection open, await the ROOT plots

            $('#processing').remove();
            $('#completed').html("database search completed: " + data.completed + "/" + va_count + ", awaiting plots...");

            //remove the progress bar
            $(".progress").remove();

            /*if (data.exists)
                fetch_json_data();
            else
                poll_status();*/

            if (loaded == true) {
                // really we should be waiting until onloaded has executed
                // await/promise?

                if (data.exists)
                    fetch_plots();
                else
                    poll_status();
            }
        }
    }
}

function onloaded() {
    loaded = true;
}

// the HR plot
/*function fetch_hr() {
    var xmlhttp = new XMLHttpRequest();
    var url = "DATA/" + uuid + "/hr.json";

    xmlhttp.onreadystatechange = function () {
        if (xmlhttp.readyState == 4 && xmlhttp.status == 200) {
            var data = xmlhttp.response;

            var plot = [
                {
                    z: data.bins,
                    type: 'heatmap'
                }
            ];

            Plotly.newPlot('hr', plot);

            var html = '<p>M<SUB>G</SUB> = phot_g_mean_mag + 5 + 5 log<SUB>10</SUB>(parallax / 1000)';

            if (where != "")
                html += ', where ' + where;

            html += '</p>';

            $('#mg').append(html);
        }
    }

    xmlhttp.open("GET", url, true);
    xmlhttp.responseType = 'json';
    xmlhttp.timeout = 0;
    xmlhttp.send();
}

// the XY plot
function fetch_xy() {
    var xmlhttp = new XMLHttpRequest();
    var url = "DATA/" + uuid + "/xy.json";

    xmlhttp.onreadystatechange = function () {
        if (xmlhttp.readyState == 4 && xmlhttp.status == 200) {
            var data = xmlhttp.response;

            var plot = [
                {
                    z: data.bins,
                    type: 'heatmap'
                }
            ];

            Plotly.newPlot('xy', plot);
        }
    }

    xmlhttp.open("GET", url, true);
    xmlhttp.responseType = 'json';
    xmlhttp.timeout = 0;
    xmlhttp.send();
}

// the RZ plot
function fetch_rz() {
    var xmlhttp = new XMLHttpRequest();
    var url = "DATA/" + uuid + "/rz.json";

    xmlhttp.onreadystatechange = function () {
        if (xmlhttp.readyState == 4 && xmlhttp.status == 200) {
            var data = xmlhttp.response;

            var plot = [
                {
                    z: data.bins,
                    type: 'heatmap'
                }
            ];

            Plotly.newPlot('rz', plot);
        }
    }

    xmlhttp.open("GET", url, true);
    xmlhttp.responseType = 'json';
    xmlhttp.timeout = 0;
    xmlhttp.send();
}

function fetch_json_data() {
    console.log("fetching json data for " + uuid);

    $('#completed').remove();

    fetch_hr();

    fetch_xy();

    fetch_rz();
}*/

function fetch_plots() {
    console.log("fetching plots for " + uuid);

    $('#completed').remove();
    //JSROOT.gStyle.fOptLogz = 1;
    JSROOT.gStyle.fPadLeftMargin = 0.15;
    JSROOT.gStyle.fPadRightMargin = 0.2;

    new JSROOT.TFile("DATA/" + uuid + "/hr.root", function (file) {
        try {
            file.ReadObject(uuid + "::HR", function (obj) {
                // revert the Y axis
                JSROOT.draw("hr", obj, "colz_ry;logz");

                var html = '<p>M<SUB>G</SUB> = phot_g_mean_mag + 5 + 5 log<SUB>10</SUB>(parallax / 1000)';

                if (where != "")
                    html += ', where ' + where;

                html += '</p>';

                $('#mg').append(html);
            });
        }
        catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#hr').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/xy.root", function (file) {
        try {
            file.ReadObject(uuid + "::XY", function (obj) {
                JSROOT.draw("xy", obj, "colz;logz");
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#xy').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/rz.root", function (file) {
        try {
            file.ReadObject(uuid + "::RZ", function (obj) {
                JSROOT.draw("rz", obj, "colz;logz");
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#rz').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/xyvr_mean.root", function (file) {
        try {
            file.ReadObject(uuid + "::XYVR_mean", function (obj) {
                JSROOT.draw("xyvr_mean", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#xyvr').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/xyvr_error.root", function (file) {
        try {
            file.ReadObject(uuid + "::XYVR_error", function (obj) {
                JSROOT.draw("xyvr_error", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#xyvr').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/xyvphi_mean.root", function (file) {
        try {
            file.ReadObject(uuid + "::XYVPhi_mean", function (obj) {
                JSROOT.draw("xyvphi_mean", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#xyvphi').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/xyvphi_error.root", function (file) {
        try {
            file.ReadObject(uuid + "::XYVPhi_error", function (obj) {
                JSROOT.draw("xyvphi_error", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#xyvphi').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/xyvz_mean.root", function (file) {
        try {
            file.ReadObject(uuid + "::XYVZ_mean", function (obj) {
                JSROOT.draw("xyvz_mean", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#xyvz').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/xyvz_error.root", function (file) {
        try {
            file.ReadObject(uuid + "::XYVZ_error", function (obj) {
                JSROOT.draw("xyvz_error", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#xyvz').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/rzvr_mean.root", function (file) {
        try {
            file.ReadObject(uuid + "::RZVR_mean", function (obj) {
                JSROOT.draw("rzvr_mean", obj);
            });
        }
        catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#rzvr').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/rzvr_error.root", function (file) {
        try {
            file.ReadObject(uuid + "::RZVR_error", function (obj) {
                JSROOT.draw("rzvr_error", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#rzvr').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/rzvphi_mean.root", function (file) {
        try {
            file.ReadObject(uuid + "::RZVPhi_mean", function (obj) {
                JSROOT.draw("rzvphi_mean", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#rzvphi').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/rzvphi_error.root", function (file) {
        try {
            file.ReadObject(uuid + "::RZVPhi_error", function (obj) {
                JSROOT.draw("rzvphi_error", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#rzvphi').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/rzvz_mean.root", function (file) {
        try {
            file.ReadObject(uuid + "::RZVZ_mean", function (obj) {
                JSROOT.draw("rzvz_mean", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#rzvz').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    new JSROOT.TFile("DATA/" + uuid + "/rzvz_error.root", function (file) {
        try {
            file.ReadObject(uuid + "::RZVZ_error", function (obj) {
                JSROOT.draw("rzvz_error", obj);
            });
        } catch (err) {
            console.log("data not found");
            $("#plots").remove();
            $('#rzvz').remove();
            $('#no-data').html("No results have been found. Please try other search criteria.");
        }
    });

    // TEST
    /*new JSROOT.TFile("DATA/" + uuid + "/rzvphi_mean.root", function (file) {
        file.ReadObject(uuid + "::RZVPhi_mean", function (obj) {
            JSROOT.draw("large", obj, "contz");
        });
    });

    new JSROOT.TFile("DATA/" + uuid + "/rzvphi_mean.root", function (file) {
        file.ReadObject(uuid + "::RZVPhi_mean", function (obj) {
            JSROOT.draw("small", obj, "contz");
        });
    });*/
}

function main() {
    session_data = document.getElementById('session-data');
    uuid = session_data.getAttribute('data-uuid');
    where = session_data.getAttribute('data-where');
}