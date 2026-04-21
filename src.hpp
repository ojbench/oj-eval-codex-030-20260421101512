#pragma once
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <utility>

typedef std::vector<std::vector<double>> IMAGE_T;

namespace nr_internal {
struct BBox {int r0,c0,r1,c1;};
static std::vector<std::vector<int>> binarize(const IMAGE_T &img) {
    int n = (int)img.size(); if(n==0) return {}; int m=(int)img[0].size();
    double sum=0, cnt=0; for(int i=0;i<n;i++) for(int j=0;j<m;j++){ sum+=img[i][j]; cnt++; }
    double mean = cnt? (sum/cnt):0.0; double thr = std::clamp(mean*0.9, 0.3, 0.7);
    std::vector<std::vector<int>> b(n, std::vector<int>(m));
    for(int i=0;i<n;i++) for(int j=0;j<m;j++) b[i][j] = (img[i][j] >= thr) ? 1 : 0;
    return b;
}

static BBox bbox(const std::vector<std::vector<int>> &b){
    int n=b.size(); if(n==0) return {0,0,0,0}; int m=b[0].size();
    int r0=n, c0=m, r1=-1, c1=-1;
    for(int i=0;i<n;i++) for(int j=0;j<m;j++) if(b[i][j]){ r0=std::min(r0,i); c0=std::min(c0,j); r1=std::max(r1,i); c1=std::max(c1,j);} 
    if(r1==-1) return {0,0,0,0};
    return {r0,c0,r1,c1};
}

static void crop(const std::vector<std::vector<int>> &b, const BBox &bb, std::vector<std::vector<int>> &out){
    if(bb.r1<bb.r0 || bb.c1<bb.c0){ out={{}}; return; }
    int h = bb.r1-bb.r0+1, w = bb.c1-bb.c0+1; out.assign(h, std::vector<int>(w,0));
    for(int i=0;i<h;i++) for(int j=0;j<w;j++) out[i][j] = b[bb.r0+i][bb.c0+j];
}

static int count_holes(const std::vector<std::vector<int>> &fg){
    int n=fg.size(); if(n==0) return 0; int m=fg[0].size();
    std::vector<std::vector<int>> bg(n+2, std::vector<int>(m+2,1));
    for(int i=0;i<n;i++) for(int j=0;j<m;j++) bg[i+1][j+1] = fg[i][j]?0:1;
    std::queue<std::pair<int,int>> q; std::vector<std::vector<int>> vis(n+2,std::vector<int>(m+2));
    auto push=[&](int i,int j){ vis[i][j]=1; q.push({i,j}); };
    for(int i=0;i<n+2;i++) for(int j=0;j<m+2;j++) if((i==0||j==0||i==n+1||j==m+1) && bg[i][j]==1){ push(i,j); }
    int di[4]={1,-1,0,0}; int dj[4]={0,0,1,-1};
    while(!q.empty()){
        auto [x,y]=q.front(); q.pop();
        for(int k=0;k<4;k++){
            int nx=x+di[k], ny=y+dj[k];
            if(nx>=0&&nx<n+2&&ny>=0&&ny<m+2 && !vis[nx][ny] && bg[nx][ny]==1) push(nx,ny);
        }
    }
    int holes=0; vis.assign(n+2,std::vector<int>(m+2));
    for(int i=1;i<=n;i++) for(int j=1;j<=m;j++) if(bg[i][j]==1 && !vis[i][j]){
        holes++; std::queue<std::pair<int,int>> qq; qq.push({i,j}); vis[i][j]=1;
        while(!qq.empty()){
            auto [x,y]=qq.front(); qq.pop();
            for(int k=0;k<4;k++){
                int nx=x+di[k], ny=y+dj[k];
                if(nx>=1&&nx<=n&&ny>=1&&ny<=m && !vis[nx][ny] && bg[nx][ny]==1) {vis[nx][ny]=1; qq.push({nx,ny});}
            }
        }
    }
    return holes;
}

static std::pair<double,double> centroid(const std::vector<std::vector<int>> &b){
    int n=b.size(); if(n==0) return {0,0}; int m=b[0].size();
    double sx=0, sy=0, s=0; for(int i=0;i<n;i++) for(int j=0;j<m;j++) if(b[i][j]){ sx+=i+0.5; sy+=j+0.5; s+=1; }
    if(s==0) return {0.5,0.5}; return {sx/s/(double)n, sy/s/(double)m};
}

static std::pair<double,double> hole_centroid(const std::vector<std::vector<int>> &fg){
    int n=fg.size(); if(n==0) return {0.5,0.5}; int m=fg[0].size();
    std::vector<std::vector<int>> bg(n, std::vector<int>(m));
    for(int i=0;i<n;i++) for(int j=0;j<m;j++) bg[i][j] = fg[i][j]?0:1;
    std::queue<std::pair<int,int>> q; std::vector<std::vector<int>> vis(n,std::vector<int>(m));
    auto push=[&](int i,int j){ vis[i][j]=1; q.push({i,j}); };
    for(int i=0;i<n;i++){ if(bg[i][0]) push(i,0); if(bg[i][m-1]) push(i,m-1);} 
    for(int j=0;j<m;j++){ if(bg[0][j]) push(0,j); if(bg[n-1][j]) push(n-1,j);} 
    int di[4]={1,-1,0,0}; int dj[4]={0,0,1,-1};
    while(!q.empty()){
        auto [x,y]=q.front(); q.pop();
        for(int k=0;k<4;k++){
            int nx=x+di[k], ny=y+dj[k];
            if(nx>=0&&nx<n&&ny>=0&&ny<m && !vis[nx][ny] && bg[nx][ny]==1) {vis[nx][ny]=1; q.push({nx,ny});}
        }
    }
    double sx=0, sy=0, s=0; for(int i=0;i<n;i++) for(int j=0;j<m;j++) if(bg[i][j]==1 && !vis[i][j]){ sx+=i+0.5; sy+=j+0.5; s+=1; }
    if(s==0) return {0.5,0.5}; return {sx/s/(double)n, sy/s/(double)m};
}

static void projections(const std::vector<std::vector<int>> &b, std::vector<int> &rows, std::vector<int> &cols){
    int n=b.size(); int m=b[0].size(); rows.assign(n,0); cols.assign(m,0);
    for(int i=0;i<n;i++) for(int j=0;j<m;j++) if(b[i][j]){ rows[i]++; cols[j]++; }
}

}

int judge(IMAGE_T &img){
    using namespace nr_internal;
    if(img.empty()||img[0].empty()) return 0;
    auto b = binarize(img);
    auto bb = bbox(b);
    std::vector<std::vector<int>> fg;
    crop(b, bb, fg);
    int h = (int)fg.size(); int w = h? (int)fg[0].size():0; if(h==0||w==0) return 0;
    int white=0; for(auto &r: fg) for(int v: r) white+=v; double density = white / (double)(h*w);
    auto cen = centroid(fg); double cy = cen.first; double cx = cen.second;
    int holes = count_holes(fg);
    if(holes>=2){ return 8; }
    if(holes==1){ auto hc = hole_centroid(fg); double hy=hc.first; 
        if(std::abs(cy-0.5) < 0.08 && std::abs(hy-0.5) < 0.12) return 0;
        if(hy < 0.5) return 9; else return 6;
    }
    double ar = (double)w / (double)h;
    if(ar < 0.5 && density < 0.55) return 1;
    std::vector<int> rows, cols; projections(fg, rows, cols);
    int top_band = 0; for(int i=0;i<std::min(h, (int)(h*0.25)); i++) top_band = std::max(top_band, rows[i]);
    int bottom_band = 0; for(int i=h-1;i>=std::max(0, h - (int)(h*0.25)); i--) bottom_band = std::max(bottom_band, rows[i]);
    int top_sum = 0; for(int i=0;i<(int)(h*0.4); i++) top_sum += rows[i];
    int bottom_sum = 0; for(int i=(int)(h*0.6); i<h; i++) bottom_sum += rows[i];
    int left_sum=0,right_sum=0; for(int j=0;j<(int)(w*0.5);j++) left_sum+=cols[j]; for(int j=(int)(w*0.5); j<w; j++) right_sum+=cols[j];
    if(top_band > (int)(0.7*w) && top_sum > (int)(1.2*bottom_sum) && cx > 0.55) return 7;
    int columns_two_segments=0;
    for(int j=0;j<w;j++){
        int seg=0; for(int i=0;i<h;i++){ if(fg[i][j] && (i==0 || !fg[i-1][j])) seg++; }
        if(seg>=2) columns_two_segments++;
    }
    if(columns_two_segments > w/3 && top_sum > 0 && rows[h/2] > rows[h*3/4]/2) return 4;
    bool has_top_bar = top_band > (int)(0.6*w);
    bool has_bottom_bar = bottom_band > (int)(0.6*w);
    int upper_left=0, upper_right=0, lower_left=0, lower_right=0;
    for(int i=0;i<h;i++) for(int j=0;j<w;j++) if(fg[i][j]){
        if(i < h/2 && j < w/2) upper_left++;
        else if(i < h/2) upper_right++;
        else if(j < w/2) lower_left++;
        else lower_right++;
    }
    if(has_top_bar && has_bottom_bar){ if(right_sum > left_sum*1.15) return 3; else return 2; }
    if(upper_left > upper_right*1.2 && lower_right > lower_left*1.2) return 5;
    if(right_sum > left_sum*1.25) return 3;
    if(left_sum > right_sum*1.25) return 2;
    if(cy < 0.45) return 7;
    return 3;
}
EOF}
