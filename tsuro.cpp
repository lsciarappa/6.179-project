#include <algorithm>
#include <array>
#include <cstring>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>

/* We use an ARRAY COORDINATE SYSTEM. This means that y is the first coord, x the second, and origin is in top left. */

#define grid_size 6
#define ascii_art_n 1
#define ascii_art_m 3

typedef std::pair<int, int> Direction;

class Tile {
private: 
  std::string s;
public:
  static int trans_way(int way){
    int perm[8] = {3, 6, 5, 0, 7, 2, 1, 4};
    return perm[way];
  }
  static Direction trans_dir(int way){
    Direction mapping[8] = {Direction(0, 1), Direction(-1, 0), Direction(-1, 0), Direction(0, -1), Direction(0, -1), \
					     Direction(1, 0), Direction(1, 0), Direction(0, 1)};
    return mapping[way];
  }
  Tile() : s("--------") {}
  Tile(int n);
  Tile(const Tile& t): s(t.s) {}
  Tile(std::string str) : s(str) {}
  std::string raw() const;
  bool valid() const;
  bool empty() const;
  int match(int) const;
  static Tile canonize(const Tile);  
};

struct Marker {
public:
  int y = (grid_size/2), x = 0, w = 3;
  bool valid() const {
    return 0 <= x && x < grid_size && 0 <= y && y < grid_size && 0 <= w && w < 8;
  }
};

class Game{
private:
  std::array<std::array<Tile,grid_size>,grid_size> grid;
  std::array<Marker,2> markers;
  bool player1turn;
  std::array<Tile,6> hand;
public:
  Game();
  Game(const Game& other);
  Tile gridAt(int x, int y) const;
  bool toMove() const;
  Tile inHand(int index) const;
  Tile inHand(int index, bool player1) const;
  Marker markerOf(bool player1) const;
  bool hasLost(bool player1p) const;
  void canonize();
  void advance();
  void play(int index, int nrotates);
  void draw();
  std::vector<Game> & children();
  double heuristic() const;
  friend std::ostream &operator<<(std::ostream &os, const Game &g);
  friend bool operator==(const Game &g1, const Game &g2);  
};

Game::Game(const Game & other){
  player1turn = other.player1turn;
  for(int i = 0; i < 6; i++){
    hand[i] = Tile(other.hand[i]);
  }
  for(int i = 0; i < 2; i++){
    markers[i].x = other.markers[i].x;
    markers[i].y = other.markers[i].y;
    markers[i].w = other.markers[i].w;
  }
  for(int i = 0; i < grid_size; i++){
    for(int j = 0; j < grid_size; j++){
      grid[i][j] = Tile(other.grid[i][j]);
    }
  }
}

Tile::Tile(int n){
  if(n>=6) throw std::domain_error("Tiles have at most ten ways"); // In the real game n = 4, so we use single-digit numbers only.
  s = "";
  for(int i = 0; i < 2 * n; i++){
    s += '*';
  }
  for(int i = 0; i < n; i++){
    int j = 0;
    while (s[j] != '*') j++;
    s[j] = (char)(i+(int)'0');
    std::random_device gen;
    std::uniform_int_distribution<int> distribution(0, 2 * n);
    while (s[j] != '*') j = distribution(gen);
    s[j] = (char)(i+(int)'0');    
  }
}

std::string Tile::raw() const{
  return s;
}

bool Tile::valid() const{
  if (s.length() % 2 != 0) return false;
  int n = s.length()/2;
  if (n != 4) return false; // In the real game n = 4, and other numbers don't make sense grid-wise
  if(s[0] != '0'){
    for (int i = 0; i < 2 * n; i++){
      if (s[i] != '-') return false;
    }
  }
  int seen[n];
  for (int i = 0; i < n; i++){
    seen[i] = 0;
  }
  int next = 0; 
  for(int i = 0; i < 2 * n; i++){
    int j = (int)(s[i]) - (int)'0';
    if (j >= n) return false;
    switch(seen[j]) {
    case 0:
      if (j != next) return false;
      next++;
      seen[j]++;
      break;
    case 1:
      seen[j]++;
      break;
    default:
      return false;
    }
  }
  for(int i = 0; i < n; i++){
    if (seen[i] != 2) return false;
  }
  return true;  
}

bool Tile::empty() const{
  return s[0] != '0';
}

int Tile::match(int way) const{
  if (way < 0 || way >= 8) throw std::domain_error("Undefined way in match");
  int i = 0;
  while (i == way || s[i] != s[way]) i++;
  return i;
}

//center: 2n + 1 = (4n+3 - 1)/2. First quart: n. 3rd quart: 3n + 2.

std::ostream &operator<<(std::ostream &os, const Tile &t){
  return (os << t.raw());
}

void pretty(std::vector<std::vector<std::string>> & out, const Tile &t, int x_offs, int y_offs){
  const int n = ascii_art_n, m = ascii_art_m;
  std::array<std::array<char, 4 * m + 3>, 4 * n + 3> art;
  for(int x = 0; x < 4*m+3; x++){
    for(int y = 0; y < 4*n+3; y++){
      art[y][x] = '0';
    }
  }
  if(t.empty()){
    for(int y = 0; y < 4*n+3; y++){
      for(int x = 0; x < 4*m+3; x++){
	out[y_offs + y][x_offs + x] = "\u001b[38;5;232m";
	out[y_offs + y][x_offs + x] += art[y][x];
      }
    }
    return;
  }
  const int						  \
    ways_x[8] = {4*m+2, 3*m+2, n, 0, 0, n, 3*m+2, 4*m+2}, \
    ways_y[8] = {n, 0, 0, n, 3*n+2, 4*n+2, 4*n+2, 3*n+2};
  for(int i = 0; i < 4; i++){
    double start_x = -1, start_y = -1, end_x, end_y;
    for(int j = 0; j < 8; j++){
      if(((int)'0') + i == ((int)t.raw()[j])){
	if(start_x == -1){
	  start_x = ways_x[j];
	  start_y = ways_y[j];
	}else{
	  end_x = ways_x[j];
	  end_y = ways_y[j];
	}
      }
    }
    for(int x = 0; x < 4*m+3; x++){
      for(int y = 0; y < 4*n+3; y++){
	if (std::pow( /* distance from point to line between points */
     (end_y-start_y) * x - (end_x-start_x)*y + end_x*start_y - start_x*end_y, 2)\
	    <= 0.3 * (std::pow(start_y-end_y,2) + std::pow(start_x-end_x,2))){
	  if((x == 0) + (y == 0) + (x == 4*m+2) + (y == 4*n+2) < 2){ //don't do corners
	    art[y][x] += std::pow(2,i);
	  }
	}
      }
    }
  }
  const std::string coloring[16] =
    {"232", "1", "3", "202", "25", "4", "2", "94",
     "255", "225", "226", "208", "32", "200", "155", "215"};
  for(int y = 0; y < 4*n+3; y++){
    for(int x = 0; x < 4*m+3; x++){
      std::string s = "\u001b[38;5;"+ coloring[(int)art[y][x] - (int)'0']+"m";
      out[y_offs + y][x_offs + x] = s + art[y][x];
    }
  }
}

Tile rotateCW(const Tile t){
  std::string s = t.raw(), out = "";
  if (t.empty()) return Tile();
  s = s.substr(2) + s.substr(0, 2);
  int perm[4] = {-1, -1, -1, -1}, count = 0;
  for(int i = 0; i < s.length(); i++){
    int n = ((int)s[i])-((int)'0');
    if (perm[n] == -1){
      perm[n] = count;
      count++;
    }
  }
  for(int i = 0; i < s.length(); i++){
    out += (char)(((int)'0')+perm[((int)s[i])-((int)'0')]);
  }
  return Tile(out);
}

Tile Tile::canonize(const Tile t){
  std::array<Tile, 4> candidates = {t, rotateCW(t), rotateCW(rotateCW(t)),
			 rotateCW(rotateCW(rotateCW(t)))};
  return *std::min_element(candidates.begin(), candidates.end());
}


std::ostream &operator<<(std::ostream &os, const Marker &m){
  return (os << "x: " << m.x << " y: " << m.y << " way: " << m.w);
}

Game::Game(){
  for(int x = 0; x < grid_size; x++){
    for(int y = 0; y < grid_size; y++){
      grid[y][x] = Tile();
    }
  }
  markers[0].x = 0; markers[0].y = (grid_size/2); markers[0].w = 3;
  markers[1].x = (grid_size-1); markers[1].y = (grid_size/2 - 1); markers[1].w = 7;
  player1turn = false; //player 0 goes first
  for(int i = 0 ; i < 6; i++){
    hand[i] = Tile(4);
  }
  this->canonize();
}

Tile Game::gridAt(int x, int y) const{
  return grid[y][x];
}

bool Game::toMove() const{
  return player1turn;
}

Tile Game::inHand(int index) const{
  return hand[index];
}

Tile Game::inHand(int index, bool player1) const{
  return hand[index + 3 * player1];
}

Marker Game::markerOf(bool player1) const{
  Marker out;
  out.x = markers[player1].x;  
  out.y = markers[player1].y;
  out.w = markers[player1].w;
  return out;
}

bool Game::hasLost(bool player1p) const{
  return !(markers[player1p].valid());
}

bool operator<(Tile & t, Tile & s){
  return t.raw() < s.raw();
}

void Game::canonize(){
  for(int i = 0; i < 6; i++){
    hand[i] = Tile::canonize(hand[i]);
  }
  std::sort(hand.begin(), hand.begin()+3);
  std::sort(hand.begin()+3, hand.end());
}

void Game::advance(){
  Marker & marker = markers[player1turn];
  while (marker.valid() && !grid[marker.y][marker.x].empty()){
    int new_inner_way = grid[marker.y][marker.x].match(marker.w);
    marker.w = Tile::trans_way(new_inner_way);
    marker.y += Tile::trans_dir(new_inner_way).first;
    marker.x += Tile::trans_dir(new_inner_way).second;
  }
  this->canonize();
}

void Game::play(int index, int nrotates){
  const int m = index + 3 * player1turn;
  for(int i = 0; i < nrotates; i++){
    hand[m] = rotateCW(hand[m]);
  }
  grid[markers[player1turn].y][markers[player1turn].x] = hand[m];
  hand[m] = Tile();
  advance();
  player1turn = !player1turn;
  advance();
}

void Game::draw(){
  for(int i = 0; i < 6; i++){
    if(hand[i].empty()) hand[i] = Tile(4);
  }
  this->canonize();
}

double Game::heuristic() const{//Positive is favorable for player 0.
  if(this->hasLost(false)) return -1 * std::numeric_limits<double>::infinity();
  if(this->hasLost(true)) return std::numeric_limits<double>::infinity();
  return ( (markers[1].x-(1.0+grid_size)/2)*(markers[1].x-(1.0+grid_size)/2)
	   + (markers[1].y-(1.0+grid_size)/2)*(markers[1].y-(1.0+grid_size)/2)
	   - (markers[0].x-(1.0+grid_size)/2)*(markers[0].x-(1.0+grid_size)/2)
	   - (markers[0].y-(1.0+grid_size)/2)*(markers[0].y-(1.0+grid_size)/2));
} // JB's idea: sum_{empty tiles} min distance from you to tile - '' opponent.

std::string pretty(const Game &g){
  std::vector<std::vector<std::string>> all, hnds;
  const int n = ascii_art_n, m = ascii_art_m, 
  width = 6 * (4 * ascii_art_m + 4)+1, height = 6 * (4 * ascii_art_n + 4)+1;
  all.resize(height);
  for(int y = 0; y < height; y++){
    all[y].resize(width);
  }
  for(int x = 0; x < width; x++){
    for(int y = 0; y < height; y++){
      all[y][x] = " ";
    }
  }

  for(int mx = 0; mx < 6; mx++){
    for(int my = 0; my < 6; my++){
      pretty(all, g.gridAt(mx, my),
	     1+ mx * (4 * ascii_art_m + 4),
	     1+ my * (4 * ascii_art_n + 4));
      
    }
    
  }
  const int marker_xs[8] = {4*m+3, 3*m+2, m, -1, -1, m, 3*m+2, 4*m+3},
    marker_ys[8] = {n, -1, -1, n, 3*n+2, 4*n+3, 4*n+3, 3*n+2};
  const char marker_syms[8] = {'<', 'v', 'v', '>', '>', '^', '^', '<'};
  const Marker m0 = g.markerOf(false), m1 = g.markerOf(true);

  for(int mx = 0; mx < 6; mx++){
    for(int my = 0; my < 6; my++){
      for(int w = 0; w < 8; w++){
	int hx = 1 + mx * (4 * ascii_art_m + 4) + marker_xs[w],
	  hy = 1 + my * (4 * n + 4) + marker_ys[w];
	if(all[hy][hx] == " " || all[hy][hx] == "\u001b[38;5;55m#"){
	  if((mx == m0.x && my == m0.y && w == m0.w) || \
	     (mx == m1.x && my == m1.y && w == m1.w)){
	    all[hy][hx] = "\u001b[38;5;255m";
	    all[hy][hx] += marker_syms[w];
	  } else {
	    all[hy][hx] = "\u001b[38;5;55m#";
	  }
	}
      }
    }
  }

  hnds.resize(4*ascii_art_n+3);
  
  for(int y = 0; y < 4*ascii_art_n+3; y++){
    hnds[y].resize(6*(4*ascii_art_m+4)+1, " ");
    hnds[y][3*(4*ascii_art_m+4)] = "\x1b" "[0m";
    hnds[y][3*(4*ascii_art_m+4)] += '|';
  }
  
  for(int h = 0; h < 6; h++){
    pretty(hnds, g.inHand(h),
	   h*(4*ascii_art_m+4)+1, 0);
  }
  
  std::string out;
  for(int y = 0; y < height; y++){
    for(int x = 0; x < width; x++){
      out += all[y][x];
    }
    out += '\n';
  }
  out += "\x1b" "[0m";
  for(int i = 0; i < 6*(4*ascii_art_m+4)+1; i++){
    out += "-";
  }
  out += '\n';

  for(int y = 0; y < 4*ascii_art_n+3; y++){
    for(int x = 0; x < 6*(4*ascii_art_m+4)+1; x++){
      out += hnds[y][x];
    }
    out += '\n';
  }
  
  out += "\x1b" "[0m";

  out += "Player ";
  out += ('0'+g.toMove());
  out += " to move.";
  return out;
}


std::ostream &operator<<(std::ostream &os, const Game &g) {
  const std::string
    upper[8] = {"         <", "        v ", " v        ", ">         ",\
		"          ", "          ", "          ", "          "},
    lower[8] = {"          ", "          ", "          ", "          ",\
		">         ", " ^        ", "        ^ ", "         <"};
  std::string outstring;
  for (int y = 0; y < grid_size; y++){
    std::string line0, line1, line2;
    for(int x = 0; x < grid_size; x++){
      line1 += ' ' + g.grid[y][x].raw() + ' ';
      if((g.markers[true].x == x && g.markers[true].y == y)){
	if((g.markers[false].x == x && g.markers[false].y == y)){
	  for(int i = 0; i < 10; i++){
	    if(upper[g.markers[true].w][i] != ' ') {
	      line0 += upper[g.markers[true].w][i];
	    } else if (upper[g.markers[false].w][i] != ' ') {
	      line0 += upper[g.markers[false].w][i];	      
	    } else {
	      line0 += ' ';
	    }
	    if(lower[g.markers[true].w][i] != ' ') {
	      line2 += lower[g.markers[true].w][i];
	    } else if (lower[g.markers[false].w][i] != ' ') {
	      line2 += lower[g.markers[false].w][i];	      
	    } else {
	      line2 += ' ';
	    }
	  }
	} else {
	  line0 += upper[g.markers[true].w];
	  line2 += lower[g.markers[true].w];	  
	}
      } else {
	if((g.markers[false].x == x && g.markers[false].y == y)){
	  line0 += upper[g.markers[false].w];
	  line2 += lower[g.markers[false].w];
	} else {
	  line0 += lower[0];
	  line2 += lower[0];
	}
      }
    }
    outstring += line0 + '\n' + line1 + '\n' + line2 + '\n';
  }

  os << outstring;
  os << "Player " << g.player1turn << " to move." << std::endl;
  os << "Player 0 at " << g.markers[0] << ", player 1 at " << g.markers[1]
     << "." << std::endl;
  os << "Player 0 hand " << g.hand[0] << ", " << g.hand[1] << ", " << g.hand[2] << "; "
     << "player 1 hand " << g.hand[3] << ", " << g.hand[4] << ", " << g.hand[5] << "." << std::endl;
  return os;
}

bool operator==(const Marker &m1, const Marker &m2){
  return m1.x == m2.x && m1.y == m2.y && m1.w == m2.w; // hi mollie!
}

bool operator==(const Tile &t1, const Tile &t2){
  return Tile::canonize(t1).raw() == Tile::canonize(t2).raw();
}

bool operator==(const Game &g1, const Game &g2){
  return g1.grid == g2.grid && g1.markers == g2.markers && \
    g1.player1turn == g2.player1turn && g1.hand == g2.hand;
}

struct gamehasher{
  size_t operator()(const Game& game) const{
    std::stringstream s;
    s << game;
    return std::hash<std::string>()(s.str());
  }
};

std::vector<Game> & Game::children(){
  std::vector<Game> * out = new std::vector<Game>;
  if (this->hasLost(false) || this->hasLost(true)){
    return *out;
  }
  for(int index = 0; index < 3; index++){
    if(!inHand(index,player1turn).empty()){
      for(int nrotates = 0; nrotates < 4; nrotates++){
	Game child = Game(*this);
	child.play(index, nrotates);
	out->push_back(child);
      }
    }
  }
  return *out;
}

double min_max(std::unordered_map<Game, double, gamehasher> & memo, Game & node, int depth){
  //  std::cout << depth;
  node.canonize();
  // if(memo.count(node)){//if value has already been computed
  //   std::cout << "already" << memo.count(node) << std::endl;
  //   return memo[node];
  // }
  if(depth > 2){//hack-y...
    memo[node] = node.heuristic();
    return node.heuristic();
  }

  
  std::vector<Game> children = node.children();
  if(children.size() == 0){
    memo[node] = node.heuristic();
    return node.heuristic();
  }
  for(auto child : children){
    min_max(memo, child,depth+1);
  }
  if(node.toMove()){ // minimize
    double out = std::numeric_limits<double>::infinity();
    for(auto child : children){
      if(memo[child] < out) out = memo[child];
    }
    memo[node] = out;
    return out;
  } else { //maximize
    double out = -1 * std::numeric_limits<double>::infinity();
    for(auto child : children){
      if(memo[child] > out) out = memo[child];
    }
    memo[node] = out;
    return out;
  }
}

int main(){
  int nhumans;
  std::cin >> nhumans;
  if(nhumans == 2){
    Game game;
    while (!game.hasLost(false) && !game.hasLost(true)){
      std::cout << pretty(game) << std::endl;
      int index, turns;
      scanf("%d %d", &index, &turns);
      game.play(index, turns);
      game.draw();
    }
    std::cout << pretty(game) << std::endl;
    if(game.hasLost(false)){
      std::cout << "Player 0 lost!" << std::endl;
    }
    if(game.hasLost(true)){
      std::cout << "Player 1 lost!" << std::endl;
    }
  } else {
    std::cout << "vs. computer" << std::endl;
    Game game;
    std::unordered_map<Game, double, gamehasher> tree = std::unordered_map<Game,double, gamehasher>();
    while (!game.hasLost(false) && !game.hasLost(true)){
      std::cout << pretty(game) << std::endl;
      min_max(tree, game,0);
      std::cout << "done" << std::endl;
      for(int index = 0; index < 3; index++){
	for(int nrotates = 0; nrotates < 4; nrotates++){
	  Game child = Game(game);
	  child.play(index, nrotates);
	  std::cout<<"tree access"<<std::endl;
	  if(tree[child] == tree[game]){
	    game.play(index, nrotates);
	    goto human_turn;
	  }
	}
      }
      std::cout << "error: no move selected" << std::endl;
    human_turn:
      if (game.hasLost(false) || game.hasLost(true)){
	goto end;
      }
      
      game.draw();
      std::cout << pretty(game) << std::endl;      
      int index, turns;
      scanf("%d %d", &index, &turns);
      game.play(index, turns);
      game.draw();
    }
  end: std::cout << pretty(game) << std::endl;
    if(game.hasLost(false)){
      std::cout << "Player 0 lost!" << std::endl;
    }
    if(game.hasLost(true)){
      std::cout << "Player 1 lost!" << std::endl;
    }
  }
  return 0;
}
