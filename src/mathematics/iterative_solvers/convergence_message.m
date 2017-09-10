function msg = convergence_message(msg, iternum, resid)

    % fprintf(repmat('\b', [1, length(msg)]));
    msg = [sprintf('iteration number: %i, relative residual: %f', iternum, resid)];
    fprintf("%s\n", msg);

end
