library_name = "libqpalm"

# Get current operating system
library_extension =
if Sys.islinux()
    "so"
elseif Sys.isapple()
    "dylib"
elseif Sys.iswindows()
    "dll"
else
    error("The platform is not supported")
end

library_filename = "$library_name.$library_extension"

library_path = joinpath(dirname(@__FILE__), "lib", library_filename)

code = "# This file is auto-generated: *DO NOT* edit, *DO NOT* commit

using Libdl

const LIBQPALM_PATH = \"$library_path\"

if Libdl.dlopen_e(LIBQPALM_PATH) == C_NULL
    error(\"Unable to load \\n\\n$library_name ($library_path)\\n\",
          \"Please make sure the library is located at the right path.\")
end
"

io = open(joinpath(dirname(@__FILE__), "deps.jl"), "w")
println(io, code)
close(io)
