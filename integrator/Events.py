def event(solution):
# Impact Definition
    breaker = True
    if solution.t_events[-1]:
        print('Impact!')
        breaker = False
    elif solution.t_events[-2]: # If a switch between apoapsis and periapsis occur, break loop
        breaker = False
        print('Periapsis Greater than Apoapsis!')
    return breaker