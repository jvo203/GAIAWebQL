function get_js_version() {
    return "JS2020-01-16.2";
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
            fetch_plots();
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

            if (loaded == true) {
                // really we should be waiting until onloaded has executed
                // await/promise?
                poll_status();
            }
        }
    }
}

function onloaded() {
    loaded = true;
}

function fetch_plots() {
    console.log("fetching plots for " + uuid);

    $('#completed').remove();
    JSROOT.gStyle.fOptLogz = 1;

    new JSROOT.TFile("DATA/" + uuid + "/hr.root", function (file) {
        file.ReadObject(uuid + "::HR", function (obj) {
            // draw only axes and revert the Y axis
            JSROOT.draw("hr", obj, "AXIS_RY");
            // overlay the actual histogram on top of the axes
            JSROOT.draw("hr", obj, "COLZ");

            var html = '<p>M<SUB>G</SUB> = phot_g_mean_mag + 5 + 5 log<SUB>10</SUB>(parallax / 1000)';

            if (where != "")
                html += ', where ' + where;

            html += '</p>';

            $('#mg').append(html);
        });
    });

    new JSROOT.TFile("DATA/" + uuid + "/xy.root", function (file) {
        file.ReadObject(uuid + "::XY", function (obj) {
            JSROOT.draw("xy", obj, "colz");
        });
    });

    new JSROOT.TFile("DATA/" + uuid + "/rz.root", function (file) {
        file.ReadObject(uuid + "::RZ", function (obj) {
            JSROOT.draw("rz", obj, "colz");
        });
    });
}

function main() {
    session_data = document.getElementById('session-data');
    uuid = session_data.getAttribute('data-uuid');
    where = session_data.getAttribute('data-where');
}