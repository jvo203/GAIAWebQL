function get_js_version() {
    return "JS2019-12-27.0";
}

function main() {
    JSROOT.gStyle.fOptLogz = 1;
    session_data = document.getElementById('session-data');
    uuid = session_data.getAttribute('data-uuid');
    where = session_data.getAttribute('data-where');

    console.log("fetching plots for " + uuid);

    new JSROOT.TFile("DATA/" + uuid + "/hr.root", function (file) {
        file.ReadObject(uuid + "::HR", function (obj) {
            JSROOT.draw("hr", obj, "colz");
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