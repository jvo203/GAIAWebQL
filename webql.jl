using Healpix
using HTTP

const halfpi=1.570796326794896619231321691639751442099
const HTMfunc_Pr=3.1415926535897932385E0/180.0

nside = 16
res = Healpix.Resolution(nside)
println(res)

function search_gaia_db(ra, dec, radius)
    θ = halfpi - dec * HTMfunc_Pr
    ϕ = ra * HTMfunc_Pr
    hpx = Healpix.ang2pixRing(res, θ, ϕ)
    println("The pixel index is $hpx")
    return "GAIA DR2 WebQL::<ra:$ra, dec:$dec, radius:$radius> --> hpx:$(hpx-1)"
end

server = "127.0.0.1"
port = 8081
htdocs = "htdocs/"

println("server: ", server, " port: ", port, " htdocs: ", htdocs)

HTTP.serve(server, port) do request::HTTP.Request
    #@show request
    #@show request.method
    #@show HTTP.header(request, "Content-Type")
    #@show HTTP.payload(request)
    println("incoming request for ", request.target)
    try               
        if request.target == "/"
            return HTTP.Response(200, read(htdocs * "index.html"))
        else
            file = htdocs * HTTP.unescapeuri(request.target[2:end])
            if isfile(file)
                HTTP.Response(200, read(file))
            else
                if occursin("GAIAWebQL.html?", request.target)
                    pos = findlast(isequal('?'), request.target)
                    params = split(request.target[pos+1:end], "&")
                    println("params:", params)
                    ra = nothing
                    dec = nothing
                    radius = nothing
                    for pair in params                        
                        if occursin("=", pair)
                            tmp = split(pair, "=")
                            key = tmp[1]
                            value = tryparse(Float64, tmp[2])
                            if key == "ra"
                                ra = value
                            end
                            if key == "dec"
                                dec = value
                            end
                            if key == "radius"
                                radius = value
                            end
                        end
                    end
                    if ra != nothing && dec != nothing && radius != nothing
                        resp = search_gaia_db(ra, dec, radius)
                        return HTTP.Response(resp)
                    else
                        return HTTP.Response(404)
                    end
                else
                    HTTP.Response(404)               
                end
            end
        end
    catch e
        return HTTP.Response(404, "Error: $e")
    end
 end

exit()