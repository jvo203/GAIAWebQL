function get_js_version() {
    return "JS2019-12-23.0";
}

/*function open_websocket_connection() {
    if ("WebSocket" in window) {
        // Let us open a web socket
        var loc = window.location, ws_uri;

        ws_uri = WS_SOCKET + loc.hostname + ':' + loc.port + "/gaiawebql/websocket/" + encodeURIComponent(uuid);

        gaia_ws = new ReconnectingWebSocket(ws_uri, null, { binaryType: 'arraybuffer' });
        gaia_ws.binaryType = 'arraybuffer';

        gaia_ws.addEventListener("open", function (evt) {
            gaia_ws.binaryType = 'arraybuffer';

            console.log("opened a websocket session");
        });

        gaia_ws.addEventListener("error", function (evt) {
            console.log('websocket conn. error:', evt);
        });

        gaia_ws.addEventListener("close", function (evt) {
            console.log("a websocket connection closed");
        });

        gaia_ws.addEventListener("message", function (evt) {
            var received_msg = evt.data;

            if (typeof evt.data === "string") {
                try {
                    var data = JSON.parse(received_msg);

                    if (data.type == "progress") {
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
                        }
                    }

                    if (data.type == "plot") {
                        //console.log(data.thread + "/" + data.total);
                        plots_received++;

                        //g_string_append_printf(gpage, "<img alt=\"search results\" src=\"data:image/png;base64,%s\" />", output) ;
                        $('#main').append('<img alt="density plot" id="' + data.thread + '" src="data:image/png;base64,' + data.density_plot + '" />');

                        $('#completed').html("awaiting plots... " + plots_received + "/" + plots_total);

                        if (data.thread == 1) {

                            var html = '<p>M<SUB>G</SUB> = phot_g_mean_mag + 5 + 5 log<SUB>10</SUB>(parallax / 1000)';

                            if (where != "")
                                html += ', where ' + where;

                            html += '</p>';

                            $('#main').append(html);
                        }

                        if (plots_received == plots_total/*data.total*/) {
    $('#processing').remove();
    $('#completed').remove();
    $(".progress").remove();
    gaia_ws.close();
}
                    }
                } catch (e) {
    console.error(e);
}
            }

if (evt.data instanceof ArrayBuffer) {
    var dv = new DataView(received_msg);
}
        });
    }
}* /

