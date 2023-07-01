import random
# This project is about a simple two-player game. In the game each player has only two choices: they can cooperate or cheat. They choose secretly, so they do not know what their opponent will do. Each player gets a score (or "payoff") as follows:

# if both players cheat, they both get 0 points

# if both players cooperate, they both get 2 points

# if one player cooperates and one player cheats, the one who cheats gets 3 points and the one who cooperates gets -1 points.

# The game is played several times and whoever has the most points at the end is the winner.

# Many different versions of this game exist, and they form a fundamental example in game theory. You can learn more about game theory in the fourth year module MATH0082.

# In this project you will investigate how well various strategies do when they play this game. You will simulate games with just two players as well as "tournaments" where many different players all play each other. Finally you will investigate the results of changing the payoffs on how well strategies perform.

# Throughout this project we will represent cooperating by the number 0 and cheating by the number 1.

cooperate = 0
cheat = 1

# A playing strategy will be represented by a function which takes a single argument and outputs either cooperate or cheat.  The argument to the function will be a list history of the choices the opponent has made so far, starting with the first choice and ending with the most recent choice. For example, if no games have been played yet then history would be the empty list [], and if three games have been played and the opponent cooperated in the first game and then cheated in the next two games, history would be equal to
# [cooperate, cheat, cheat]

def always_cheat(history):
    return cheat

# Write functions implementing the following strategies.

# Always cooperate.
# ("tit for tat") On the first game, cooperate. After that do whatever the opponent did in the previous game.
# ("random choice") In each game, cheat with probably 1/2, otherwise cooperate.
# ("grudge") If your opponent has ever cheated, cheat. Otherwise cooperate.

# For part 3, it may be helpful to import random and investigate random.choice or use to use the function random.random() which generates a uniform random number between 0 and 1.

def always_cooperate(history):
    return cooperate

def tit_for_tat(history):
    """
    Input: list containing history of choices from games played
    Output: If first game, return cooperate. If not, do what opponent did last game
    
    """
    if history == []:
        return cooperate
    else:
        return history[-1]

def random_choice(history):
    """
    Input: List containing history of choices from games played
    Output: 1/2 probability cheat, 1/2 probability cooperate

    """
    return random.choice([cheat, cooperate])

def grudge(history):
    """
    Input: List containing history of choices from games played
    Output: If list contains cheat, return cheat, ohterwise return cooperate
    
    """
    if cheat in history:
        return cheat
    else:
        return cooperate

# Write a function play_n_games(player1, player2, n). The arguments player1 and player2 should be strategy functions like the ones you wrote in Exercise 1.  The function should run n games between the two players, and should return a tuple (score1, score2) where score1 is player1's score after n games and score2 is player2's score after n games.
# You will need to keep lists history1 and history2 of the moves made by player1 and player2 to supply the arguments to the player functions, and update these lists each time a game is played.
# Print the values of:

# play_n_games(always_cooperate, tit_for_tat, 10)
# play_n_games(always_cheat, tit_for_tat, 10)
# play_n_games(grudge, tit_for_tat, 10)
# play_n_games(grudge, always_cheat, 10)

# The output of play_n_games(always_cheat, always_cooperate, 10) should be (30, -10) - if your function doesn't give this result, something is wrong.

def play_n_games(player1, player2, n):
    """
    Input: A game strategy for each player and number of times they play
    Output: A tuple showing the final score of both players
    
    """
    history1 = []     # we want to keep adding to the history of choices
    history2 = []
    score1 = 0     # we want a running total for both scores
    score2 = 0
    
    for i in range(n):    # n games are played
        a = player1(history2)
        b = player2(history1)
        score1 = score1 + 2 + a - (3 * b)     # adds 2 for a = b = 0, adds 3 if a = 1 and b = 0, subtracts 1 if a = 0 and b = 1, stays the same for a = b = 1
        score2 = score2 + 2 + b - (3 * a)     # by symmetry
        history1.append(a)     # adds the choice of cheat or cooperate to updated history list
        history2.append(b)
    return (score1, score2)

print(play_n_games(always_cooperate, tit_for_tat, 10))
print(play_n_games(always_cheat, tit_for_tat, 10))
print(play_n_games(grudge, tit_for_tat, 10))
print(play_n_games(grudge, always_cheat, 10))

# Write a function tournament(player_list).  The argument player_list will be a list of strategy functions [player1, ..., playerN].  Every player should play every other player exactly ten times - that is, you need to call play_n_games once with n=10 for each pair of distinct players from player_list.  The function should return a list [score1, ..., scoreN] of the total scores of each player.
# Print the output of

# tournament([always_cheat, always_cheat, always_cheat, always_cheat, tit_for_tat, tit_for_tat, tit_for_tat])
# tournament([always_cheat, always_cooperate, tit_for_tat, tit_for_tat])
# tournament([grudge, grudge, grudge, random_choice, always_cheat, always_cheat, tit_for_tat, tit_for_tat, tit_for_tat])

def tournament(player_list):

    totalscore = []    # later need to append to this list the scores of the n players
    for i in range(len(player_list)):
        c = player_list.copy()    # copy of list created so that we can delete elements from the copy and modify the list
        del(c[i])
        score = 0    # score for each player starts at 0
        for j in c:
            score = score + (play_n_games(player_list[i], j, 10)[0])    # the running total, no repeat games
        totalscore.append(score)    # adds the individual player score to the list
    return totalscore

print(tournament([always_cheat, always_cheat, always_cheat, always_cheat, tit_for_tat, tit_for_tat, tit_for_tat]))
print(tournament([always_cheat, always_cooperate, tit_for_tat, tit_for_tat]))
print(tournament([grudge, grudge, grudge, random_choice, always_cheat, always_cheat, tit_for_tat, tit_for_tat, tit_for_tat]))

# We will now change the payoffs so that

# if both players cheat, they both get 0 points
# if both players cooperate, they both get 1 point
# if one player cooperates and one player cheats, the one who cheats gets 3 points and the one who cooperates gets 2 points.

# This changes the game in an important way.  In the old game, no matter what player 1 does, player 2 gets a better score by cheating.  In the new game if player 1 cheats then player 2 does better by cooperating (player 2 gets 2 points by cooperating instead of 0 if they cheat), but if player 1 cooperates then player 2 does better by cheating (player 2 gets 3 points if the cheat instead of 1 if they cooperate).
# With the new payoffs, do the following:

# compute the results of 50000 games between a random_choice player 1 and a player 2 who cheats with probability p=0/10.
# compute the results of 50000 games between a random_choice player 1 and a player 2 who cheats with probability p=1/10.
# compute the results of 50000 games between a random_choice player 1 and a player 2 who cheats with probability p=2/10.
# ...
# compute the results of 50000 games between a random_choice player 1 and a player 2 who cheats with probability p=10/10.

# You should see that player 2 gets approximately the same score no matter what the value of p - this situation is called a Nash equilibrium.
# It might be helpful to use the function random.random(), which produces a uniformly distributed random number between 0 and 1.  This means the probability that random.random() is less than p equals p.

def play_n_games_ver2(prob_cheat):
    
    score1 = 0
    score2 = 0
    for x in range(50000):
        f = random.choice([cheat, cooperate])    # always random choice for player 1
        if random.random() <= prob_cheat:    # probability that the random number between 0 and 1 is less than p equals p
            g = cheat
        else:
            g = cooperate
        score1 = score1 + 1 - f + (3 * f - g) * (f - g)    # adds 1 for f = g = 0, adds 3 if f = 1 and g = 0, adds 2 if f = 0 and g = 1, stays the same for f = g = 1
        score2 = score2 + 1 - g + (3 * g - f) * (g - f)    # by symmetry
    return (score1, score2)

for y in range(11):
    prob = y/10    # multiples of 1/10 between 0 and 1 are the probabilities we use as inputs for each loop
    print(play_n_games_ver2(prob))

