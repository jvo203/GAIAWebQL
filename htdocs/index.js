function reset_gaia() {
    document.getElementById("xmin").value = "";
    document.getElementById("xmax").value = "";

    document.getElementById("ymin").value = "";
    document.getElementById("ymax").value = "";

    document.getElementById("zmin").value = "";
    document.getElementById("zmax").value = "";

    document.getElementById("rmin").value = "";
    document.getElementById("rmax").value = "";

    document.getElementById("phimin").value = "";
    document.getElementById("phimax").value = "";

    document.getElementById("mgmin").value = "";
    document.getElementById("mgmax").value = "";

    document.getElementById("parallax_over_error_min").value = "";

    document.getElementById("where").value = "";
}

function view_gaia() {
    var xmin = document.getElementById("xmin").value;
    var xmax = document.getElementById("xmax").value;

    var ymin = document.getElementById("ymin").value;
    var ymax = document.getElementById("ymax").value;

    var zmin = document.getElementById("zmin").value;
    var zmax = document.getElementById("zmax").value;

    var rmin = document.getElementById("rmin").value;
    var rmax = document.getElementById("rmax").value;

    var phimin = document.getElementById("phimin").value;
    var phimax = document.getElementById("phimax").value;

    var mgmin = document.getElementById("mgmin").value;
    var mgmax = document.getElementById("mgmax").value;

    var parallax_over_error = document.getElementById("parallax_over_error_min").value;

    var where = document.getElementById("where").value.trim();

    var params = "";

    if (xmin != "")
        params += "&xmin=" + encodeURIComponent(xmin);

    if (xmax != "")
        params += "&xmax=" + encodeURIComponent(xmax);

    if (ymin != "")
        params += "&ymin=" + encodeURIComponent(ymin);

    if (ymax != "")
        params += "&ymax=" + encodeURIComponent(ymax);

    if (zmin != "")
        params += "&zmin=" + encodeURIComponent(zmin);

    if (zmax != "")
        params += "&zmax=" + encodeURIComponent(zmax);

    if (rmin != "")
        params += "&rmin=" + encodeURIComponent(rmin);

    if (rmax != "")
        params += "&rmax=" + encodeURIComponent(rmax);

    if (phimin != "")
        params += "&phimin=" + encodeURIComponent(phimin);

    if (phimax != "")
        params += "&phimax=" + encodeURIComponent(phimax);

    if (mgmin != "")
        params += "&mgmin=" + encodeURIComponent(mgmin);

    if (mgmax != "")
        params += "&mgmax=" + encodeURIComponent(mgmax);

    if (parallax_over_error != "")
        params += "&parallax_over_error=" + encodeURIComponent(parallax_over_error);

    if (where != "")
        params += "&where=" + encodeURIComponent(where);

    if (params != "") {
        params = params.substr(1);
        var url = "/gaiawebql/GAIAWebQL.html?" + params;

        window.location.href = url;
    }

}