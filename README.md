[English version](README.en.md)

# dd_seq

蛋白質の既知の構造すべて——UniProtアクセッション番号に紐付く全PDBエントリ
（X線でもEMでも、オリゴマー状態や断片は問わない）に加えAlphaFold DBの予測
構造——を、UniProt正準配列に対する配列カバレッジと、アクティブサイト基準の
構造（RMSD）アラインメントという2軸で比較する。特定の標的に依存しない
再利用可能なパッケージとして設計されている（`dd_prep`/`dd_af`/`dd_viewer`
などと同じ方針——以下の例はすべてヒトCDK1、UniProt `P06493` を用いるが、
どのアクセッション番号でも動作する）。`dd_prep`（構造ダウンロード、HETATM
分類）と`dd_af`（fpocketベースのポケット検出）をそれぞれ再実装せず直接
再利用する。

- **Fetch（`dd_seq-fetch`）**: `list_pdb_ids_for_uniprot`（RCSB Search
  API）が、そのアクセッション番号に紐付く全PDBエントリを見つけ、各エントリ
  は`dd_prep.fetch.download_pdb`経由で、AlphaFold DBモデルは
  `dd_prep.fetch.download_afdb`経由で、正準配列はUniProt REST API経由で
  それぞれダウンロードされる。同じ`-o`ディレクトリに対して再実行すると、
  既にディスク上にあるもの（canonical.fasta、各PDBエントリ、AlphaFold
  モデル）は再ダウンロードせずスキップし、それぞれについて`already
  downloaded, skipping`と表示する——`list_pdb_ids_for_uniprot`自体は毎回
  新たに問い合わせるため、新しく公開されたエントリは他をすべて再取得する
  ことなく再実行時に拾われる。ごく最近公開された一部のエントリには
  legacyの`.pdb`ファイルがまだ存在しない（mmCIFのみ）——これらは致命的
  エラーにはせずスキップし、`manifest.json`の`"skipped"`リストに記録する。
- **Align（`dd_seq-align`）**: 取得済みの全構造について、各鎖の配列を抽出
  し（`sequence.py`、Biopython）、正準配列に対してglocalアラインメント
  （両端ギャップフリー）を行う——入力はすべて*同一*蛋白質の断片／アイソ
  フォームであり、乖離したホモログ集合ではないため、本格的なMSAツールは
  不要で、単一の参照配列で十分である。目的の蛋白質に実際に対応する鎖は
  識別性によって選ばれる（`pick_target_chain`、生のカバレッジではなく
  一致残基数でランク付け——これは、CAK-CDK1-サイクリンB1アセンブリ中の
  CDK7のような相同なパートナー鎖を誤って選ばないために必要である。この鎖
  は真の対象鎖より*高い*カバレッジを示すことがあるが、その大半はミス
  マッチである）。

  続いて、1つの「サイトソース」構造上でアクティブサイトを一度だけ定義し
  （`activesite.py`、2モード——`--site-mode ligand`: 自動選択された結合
  リガンド近傍の残基、`--site-mode pocket`: `dd_af.pocket`経由のfpocketの
  トップランクのdruggableポケット）、UniProt正準位置を経由して往復させる
  （`map_site_to_structure`）ことで他の各構造自身の残基番号に変換する——
  これにより、番号付けや鎖構成が全く異なる構造間でもサイトを比較可能に
  する。その後、全構造をPyMOL（`structalign.py`）経由で1つの参照構造に
  重ね合わせる（既定: AlphaFoldモデル。常に全長かつリガンド無しである
  ため）: `ligand`/`pocket`モードでは既知のサイト残基対応に対する
  `cmd.pair_fit`を、`--site-mode none`では`cmd.cealign`（トポロジーに
  依存しないCE構造アラインメント、残基対応不要）を用いる。サイトを全く
  解決できない構造（例: 蛋白質の折り畳みドメインではなく無関係な断片の
  周りで結晶化した共結晶）は、バッチ全体を中断せず、記録された理由と共に
  スキップされる。
- **Run（`dd_seq-run`）**: fetchとalignを1ステップで実行する。
- **App（`streamlit run app.py -- --report-dir DIR`）**: 3タブ構成——
  Overview（構造ごとのmethod/resolution/coverage/RMSD表——テーブル自体の
  内部スクロールバーが出ないよう、行数に合わせて`height`を計算して渡し、
  全行を常に表示する。ページ自体は必要に応じて普通にスクロールする）、
  Sequence
  coverage（正準位置にまたがる構造ごとのmatch/mismatch/not-resolved
  トラック）、Structure overlay（py3Dmol、各構造の対象鎖を重ね合わせて
  色分けし、アクティブサイトをハイライト、リガンド表示は任意）。表示する
  構造は、ドロップダウン式のmultiselectではなく構造ごとのチェックボックス
  （初期値は全てオン）で選ぶ——「Select all」/「Deselect all」ボタンも
  添えてあり、十数構造まとめての重ね合わせと、少数だけに絞った比較との間を
  素早く行き来できる（全オフにしてから見たい構造だけチェックし直す方が、
  多数のpillが並ぶmultiselectから1つずつ外すより速い）。3D
  シーンはプレーンな`st.components.v1.html`呼び出しではなく`dd_viewer`の
  ダブルバッファリング`view3d`コンポーネント経由で埋め込まれるため、
  ウィジェット操作のたびに既定のフィット表示へ戻るのではなく、カメラ位置
  （回転・ズーム）が保持される。「Reset view」ボタンで既定のフィット表示に
  戻せる。「Highlight active-site residues」チェックボックスは、レポートが
  `--site-mode none`で生成されている場合（ハイライトすべきサイトが無い
  場合）は無効化される。「Show only active-site surroundings」チェック
  ボックス（同じく`--site-mode none`では無効化）は、アクティブサイトを
  ハイライトするだけでなく、各構造の全体cartoonを描画せずアクティブサイト
  から指定半径（「Pocket radius」スライダー、既定8Å）以内の残基だけを
  細いスティック表示に切り替える（半径0.08——3Dmol.jsの純粋な`line`
  表示は大半のブラウザでGL_LINEの太さが1pxに固定され調整できず細すぎた
  ため、太さを制御できる細いシリンダー=スティックに変更した）——3Dmol.js
  の`within`+`byres`セレクタでサイト残基周辺の原子を残基単位で抽出し、
  それ以外は完全に非表示にする。この蛋白ワイヤーは構造ごとの色を白方向に
  薄めた色で描画し、リガンド（表示している場合）は薄めない元の色のまま
  描くことで、リガンドが背景の蛋白より視覚的に際立つようにしている。
  さらに蛋白ワイヤー・リガンド・ハイライトされたアクティブサイト残基の
  いずれも、炭素だけを（薄めた構造色／通常の構造色／`SITE_COLOR`の）
  色に、それ以外の元素（O/N/S/Pやハロゲンなど）はRasMol系の標準元素色に
  塗り分けるカラースキームを使う——3Dmol.js組み込みの`"*Carbon"`スキーム
  （例: `"yellowCarbon"`）と同じ仕組みだが、named CSS色に限定されず
  どんな16進色にも使えるよう自前で実装している。ハイライトされた残基は
  `focus_on_site`の設定に関わらず常にこの元素色分けで描画される。

## インストール

Biopython、pandas、numpy、PyMOL（`pymol2`、GUIではなくライブラリとして
import可能。conda-forgeのパッケージ名は`pymol-open-source`だが、
`pip list`/`pip show`から見えるディストリビューション名は`pymol`のため、
`pyproject.toml`では`pymol`として明記済み）、`fpocket` CLI（pipでは
入らないためconda-forge必須）、`dd_prep`/`dd_af`パッケージが必要。`dd`
conda envには既にすべて揃っている:

```bash
cd dd_prep && pip install -e . && cd ..   # if not already installed
cd dd_af && pip install -e . && cd ..     # if not already installed
cd dd_viewer && pip install -e . && cd .. # if not already installed ([app] extra のみで使用、後述)
cd dd_seq && pip install -e ".[app]"      # [app] adds streamlit/py3Dmol/matplotlib/dd_viewer
```

これにより3つのコンソールコマンド、`dd_seq-fetch`、`dd_seq-align`、
`dd_seq-run`がインストールされる。

## 使い方

```bash
dd_seq-run P06493 -o data --site-mode ligand
streamlit run app.py -- --report-dir data
```

`--site-mode`（既定`ligand`）: `ligand`（結合リガンド近傍の残基でフィット）、
`pocket`（fpocketが自動検出したdruggableポケットでフィット、apo構造や
AlphaFoldモデルにも使える）、`none`（アクティブサイトによる制限無し、
全鎖CEアラインメント）。`--reference`/`--site-source`は上記の既定値を
上書きする。`--ligand-cutoff`/`--pocket-rank`はサイト検出を調整する。

3つのコマンドはいずれも、完了したアイテムごとに1行ずつその場で出力する
（構造ごとのfetch/skip、構造ごとの配列アラインメント結果、構造ごとの
構造フィット結果またはスキップ理由）——`--no-progress`を渡すとこれを抑制
し、最終サマリー表のみを出力する。

## 設計ノート

- **本格的なMSAツールを使わない理由**: mafft/clustaloは`mpro` envに無く、
  そもそもここでは適切なツールでもない——ここでの全構造は同一の蛋白質で
  あるため、UniProt正準配列に対する参照基準のペアワイズglocalアラインメント
  （Biopythonの`PairwiseAligner`、BLOSUM62、両端ギャップフリー）の方が、
  汎用的な多重配列アラインメントよりも直接的に有用な結果（全構造にまたがる
  正準位置ごとのカバレッジ／ミスマッチ表）を一度に与える。
- **共通座標系としてのUniProt正準位置**: どの構造も独自のオーサー残基番号
  （オフセット、挿入コード、欠損密度によるギャップ）を持つ。これらの
  番号付けをペアごとに調整しようとするのではなく、すべて（アクティブ
  サイト残基、カバレッジトラック）をUniProt正準位置で表現し、使用する
  時点でのみ各構造自身の番号に変換する（`ChainAlignment.
  resseq_for_canonical`/`canonical_for_resseq`）。
- **サイトモードのフィッティングで`align`/`cealign`より`pair_fit`を使う**:
  サイト残基の対応は（両側とも同じ正準位置であるため）既に正確にわかって
  いるので、独自の内部構造／配列再マッチングを行う`cealign`/`align`では
  なく、与えられた原子対に対する直接的なKabsch重ね合わせである
  `cmd.pair_fit`を使う——PyMOLが誤った残基を黙ってペアにしてしまうリスク
  が無い。
- **3Dタブで`dd_viewer`の`view3d`コンポーネントを再利用する**: プレーンな
  `st.components.v1.html(view._make_html())`呼び出しは、Streamlitの
  再実行のたびに（実際にシーンが変わる操作かどうかに関わらず）iframeの
  中身全体を置き換えてしまうため、毎回カメラが既定のフィット位置へ
  リセットされる。`dd_viewer`は自身のpy3Dmol埋め込みで既にこの問題に
  対処済みで、小さな静的・ダブルバッファリングのStreamlitコンポーネント
  （`dd_viewer.component.view3d` + `dd_viewer.scene.
  html_with_camera_events`）を使っている——このコンポーネント自身のJSは
  （短命なシーン側iframeとは違い）再実行をまたいで生き続け、新しいシーンを
  表示する前に直前のカメラ位置を再適用する。同じ仕組みを`dd_seq`側で
  再実装するのではなく直接再利用しているのが、`dd_viewer`が`[app]`
  extraの依存関係になっている理由である。

## 既知の制約

- 目的蛋白質が短い無関係なペプチド断片のみを提供する共結晶構造
  （パートナー蛋白質自身のサイトに結合した折り畳みドメインではない場合）
  は、アクティブサイト領域を意味のあるレベルでカバーする鎖を持たない——
  `dd_seq-align`は誤ってフィットさせるのではなく、これらを正しくスキップ
  する（構造ごとの`report.json`の`"align_error"`を参照）。
- `site_from_pocket`/fpocketは、妥当なポケットを検出するために単一の
  孤立した鎖を必要とする。他の全鎖を除去した対象鎖に対して実行される
  ため、鎖間にのみ存在するポケット（例: 蛋白質間界面にのみ存在する溝）は
  この方法では見つからない。
