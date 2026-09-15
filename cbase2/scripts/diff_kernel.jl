# diff_kernel.jl — drift check for the cbase2 model kernel
#
# cbase2/src/core/ holds trimmed copies of the parent repo's working src/
# files (hybrid decision, 2026-09-15). This script compares every copied
# file against its parent-side counterpart byte-for-byte (SHA-256) and
# reports SAME / DIFF / MISSING so kernel drift is flagged, not silent.
#
# Run from anywhere:  julia cbase2/scripts/diff_kernel.jl
# Note: cbase2 divergence is legitimate once Stage 1 edits land — the
# report exists to make drift explicit, not to forbid it.

using SHA

const PARENT = normpath(joinpath(@__DIR__, "..", ".."))
const CORE   = normpath(joinpath(@__DIR__, "..", "src", "core"))
const SRC    = normpath(joinpath(@__DIR__, "..", "src"))

const CORE_FILES = [
    "interface.jl",
    "solution.jl",
    "ces.jl",
    "mobile_labor.jl",
    "leontief.jl",
    "util.jl",
    "variance_decomposition.jl",
]

function sha256_file(path::AbstractString)
    ctx = SHA.SHA2_256_CTX()
    open(path, "r") do io
        while !eof(io)
            SHA.update!(ctx, read(io, 65536))
        end
    end
    return bytes2hex(SHA.digest!(ctx))
end

clean = true
for f in CORE_FILES
    core_path   = joinpath(CORE, f)
    parent_path = joinpath(PARENT, "src", f)
    if !isfile(core_path)
        println("MISSING   core/$f  (not yet copied from src/)")
        global clean = false
    elseif !isfile(parent_path)
        println("PARENT-GONE  core/$f  (parent src/$f no longer exists)")
        global clean = false
    else
        h1, h2 = sha256_file(core_path), sha256_file(parent_path)
        println(h1 == h2 ? "SAME      core/$f" : "DIFF      core/$f  (diverged from ../src/$f)")
    end
end

# cbase2-specific files: existence check only, never compared to the parent
for f in ["financing.jl", "closures.jl", "validation.jl"]
    p = joinpath(SRC, f)
    println(isfile(p) ? "LOCAL     src/$f  (cbase2-only, no parent counterpart)" :
                         "PENDING   src/$f  (Stage 1: not yet written)")
end

exit(clean ? 0 : 1)
