using LibGit2

function audit_production_checkout(arguments::Vector{String})
    if length(arguments) != 1
        println(stderr, "usage: audit_production_checkout.jl REPOSITORY")
        return 2
    end

    repository_path = try
        realpath(only(arguments))
    catch error
        println(stderr, "checkout audit could not resolve the repository: ", error)
        return 4
    end

    repository = try
        LibGit2.GitRepo(repository_path)
    catch error
        println(stderr, "checkout audit could not open the repository: ", error)
        return 4
    end

    head_commit = ""
    origin_commit = ""
    tracked_dirty = true
    staged_dirty = true
    try
        LibGit2.need_update(repository)
        head_commit = LibGit2.head(repository_path)
        origin_commit = string(
            LibGit2.revparseid(repository, "refs/remotes/origin/DMRG")
        )
        tracked_dirty = LibGit2.isdirty(repository)
        staged_dirty = LibGit2.isdirty(repository; cached=true)
    catch error
        println(stderr, "checkout audit failed closed: ", error)
        return 4
    finally
        close(repository)
    end

    println("head_commit=", head_commit)
    println("origin_dmrg=", origin_commit)
    println("tracked_dirty=", tracked_dirty)
    println("staged_dirty=", staged_dirty)

    failed = false
    if tracked_dirty
        println(stderr, "tracked worktree changes block Fig. 2 production")
        failed = true
    end
    if staged_dirty
        println(stderr, "staged index changes block Fig. 2 production")
        failed = true
    end
    if head_commit != origin_commit
        println(stderr, "W003 HEAD must equal origin/DMRG before Fig. 2 production")
        failed = true
    end
    return failed ? 4 : 0
end

exit(audit_production_checkout(ARGS))
